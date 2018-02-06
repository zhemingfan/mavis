from copy import copy
import math
import re
import subprocess

import pysam

from . import cigar as _cigar
from .cigar import EVENT_STATES, QUERY_ALIGNED_STATES, REFERENCE_ALIGNED_STATES, convert_cigar_to_string
from ..constants import CIGAR, DNA_ALPHABET, ORIENT, READ_PAIR_TYPE, STRAND, SVTYPE, NA_MAPPING_QUALITY
from ..interval import Interval


class SamRead(pysam.AlignedSegment):
    """
    Subclass to extend the pysam.AlignedSegment class adding some utility methods and convenient representations

    Allows next_reference_name and reference_name to be set directly so that is does not depend on a bam header
    """

    def __init__(self, reference_name=None, next_reference_name=None, alignment_score=None, **kwargs):
        pysam.AlignedSegment.__init__(self)
        self._reference_name = reference_name
        self._next_reference_name = next_reference_name
        self.alignment_score = alignment_score
        self.mapping_quality = NA_MAPPING_QUALITY
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    @property
    def alignment_id(self):
        return '{}:{}[{}]{}'.format(self.reference_name, self.reference_start, self.query_name, convert_cigar_to_string(self.cigar))

    def __repr__(self):
        return '{}({}:{}-{}, {}, {}...)'.format(
            self.__class__.__name__, self.reference_name, self.reference_start, self.reference_end,
            convert_cigar_to_string(self.cigar), self.query_sequence[:10]
        )

    @property
    def query_length(self):
        return len(self.query_sequence) + sum([f for v, f in self.cigar if v == CIGAR.H] + [0])

    @classmethod
    def copy(cls, pysamread):
        return cls.copy_onto(pysamread)

    @classmethod
    def copy_onto(cls, pysamread, copyread=None):
        cp = cls() if copyread is None else copyread
        cp._reference_name = pysamread.reference_name
        cp.query_sequence = pysamread.query_sequence
        cp.reference_start = pysamread.reference_start
        cp.reference_id = pysamread.reference_id
        cp.cigar = pysamread.cigar[:]
        cp.query_name = pysamread.query_name
        cp.mapping_quality = pysamread.mapping_quality
        cp.set_tags(pysamread.get_tags())
        cp.flag = pysamread.flag
        if pysamread.is_paired:
            cp.next_reference_id = pysamread.next_reference_id
            cp.next_reference_start = pysamread.next_reference_start
            cp._next_reference_name = pysamread.next_reference_name
        try:
            cp.alignment_score = pysamread.alignment_score
        except AttributeError:
            pass
        return cp

    def __copy__(self):
        return self.__class__.copy(self)

    @property
    def reference_name(self):
        return self._reference_name

    @property
    def next_reference_name(self):
        return self._next_reference_name

    def deletion_sequences(self, reference_genome):
        """returns the reference sequences for all deletions"""
        rpos = self.reference_start
        result = []
        for state, freq in self.cigar:
            if state in REFERENCE_ALIGNED_STATES:
                if state not in QUERY_ALIGNED_STATES:
                    result.append(reference_genome[self.reference_name].seq[rpos:rpos + freq])
                rpos += freq
        return result

    def insertion_sequences(self):
        """returns the inserted sequence for all insertions"""
        qpos = 0
        result = []
        for state, freq in self.cigar:
            if state in QUERY_ALIGNED_STATES:
                if state not in REFERENCE_ALIGNED_STATES:
                    result.append(self.query_sequence[qpos:qpos + freq])
                qpos += freq
        return result


def pileup(reads, filter_func=None):
    """
    For a given set of reads generate a pileup of all reads (exlcuding those for which the filter_func returns True)

    Args:
        reads (iterable of pysam.AlignedSegment): reads to pileup
        filter_func (callable): function which takes in a  read and returns True if it should be ignored and False otherwise

    Returns:
        iterable of tuple of int and int: tuples of genomic position and read count at that position

    Note:
        returns positions using 1-based indexing
    """
    hist = {}  # genome position => frequency count
    for read in reads:
        if filter_func and filter_func(read):
            continue
        for pos in read.get_reference_positions():
            hist[pos + 1] = hist.get(pos + 1, 0) + 1
    return sorted(hist.items())


def map_ref_range_to_query_range(read, ref_range):
    """
    Args:
        ref_range (Interval): 1-based inclusive
        read (pysam.AlignedSegment): read used for the mapping
    Returns:
        Interval: 1-based inclusive range
    """
    rpos = read.reference_start
    qpos = 0
    qstart = None
    qend = None
    for state, value in read.cigar:
        for i in range(0, value):
            if state in QUERY_ALIGNED_STATES:
                qpos += 1
            if state in REFERENCE_ALIGNED_STATES:
                rpos += 1
            if qstart is None and ref_range.start <= rpos:
                qstart = qpos
            if ref_range.end >= rpos:
                qend = qpos
    if qstart is None or qend is None:
        raise ValueError('reference range is not mapped by input read', ref_range, read.reference_start, read.cigar)
    return Interval(qstart, qend)


def get_samtools_version():
    """
    executes a subprocess to try and run samtools and parse the version number from the output

    Example:
        >>> get_samtools_version()
        (1, 2, 1)
    """
    proc = subprocess.getoutput(['samtools'])
    for line in proc.split('\n'):
        match = re.search(r'Version: (?P<major>\d+)(\.(?P<mid>\d+)(\.(?P<minor>\d+))?)?', line)
        if match:
            major = int(match.group('major'))
            mid = int(match.group('mid')) if match.group('mid') else 0
            minor = int(match.group('minor')) if match.group('minor') else 0
            return major, mid, minor
    raise ValueError('unable to parse samtools version number')


def samtools_v0_sort(input_bam, output_bam):
    prefix = re.sub('\.bam$', '', output_bam)
    return 'samtools sort {} {}'.format(input_bam, prefix)


def samtools_v1_sort(input_bam, output_bam):
    return 'samtools sort {} -o {}'.format(input_bam, output_bam)


def breakpoint_pos(read, orient=ORIENT.NS):
    """
    assumes the breakpoint is the position following softclipping on the side with more
    softclipping (unless and orientation has been specified)

    Args:
        read (:class:`~pysam.AlignedSegment`): the read object
        orient (ORIENT): the orientation

    Returns:
        int: the position of the breakpoint in the input read
    """
    typ, freq = read.cigar[0]
    end_typ, end_freq = read.cigar[-1]
    ORIENT.enforce(orient)

    if typ != CIGAR.S and end_typ != CIGAR.S:
        raise AttributeError('cannot compute breakpoint for a read without soft-clipping', read.cigar)

    if orient == ORIENT.NS:
        if (typ == CIGAR.S and end_typ == CIGAR.S and freq > end_freq) \
                or typ == CIGAR.S and end_typ != CIGAR.S:
            orient = ORIENT.RIGHT
            # soft clipped to the left
        else:
            # soft clipped to the right
            orient = ORIENT.LEFT

    if orient == ORIENT.RIGHT:
        if typ != CIGAR.S:
            raise AttributeError('soft clipping doesn\'t support input orientation for a breakpoint', repr(orient), read.cigar, read.get_tags())
        return read.reference_start
    else:
        if end_typ != CIGAR.S:
            raise AttributeError('soft clipping doesn\'t support input orientation for a breakpoint', orient, read.cigar, read.get_tags())
        return read.reference_end - 1


def calculate_alignment_score(read, consec_bonus=1):
    """
    calculates a score for comparing alignments

    Args:
        read (pysam.AlignedSegment): the input read

    Returns:
        float: the score
    """
    score = 0
    qlen = read.reference_end - read.reference_start
    max_score = qlen + (qlen - 1) * consec_bonus
    for c, v in read.cigar:
        if c == CIGAR.M:
            raise ValueError('cannot calculate the alignment score if mismatch v match has not been specified')
        elif c == CIGAR.EQ:
            score += v + (v - 1) * consec_bonus
    return score / max_score


def nsb_align(
        ref, seq,
        weight_of_score=0.5,
        min_overlap_percent=1,
        min_match=0,
        min_consecutive_match=1,
        scoring_function=calculate_alignment_score):
    """
    given some reference string and a smaller sequence string computes the best non-space-breaking alignment
    i.e. an alignment that does not allow for indels (straight-match). Positions in the aligned segments are
    given relative to the length of the reference sequence (1-based)

    Args:
        ref (str): the reference sequence
        seq (str): the sequence being aligned
        weight_of_score (float): when scoring alignments this determines the amount
            of weight to place on the cigar match. Should be a number between 0 and 1
        min_overlap_percent (float): the minimum amount of overlap of the input sequence to the reference
            should be a number between 0 and 1
        min_match (float): the minimum number of matches compared to total
        scoring_function (callable): any function that will take a read as input and return a float
          used in comparing alignments to choose the best alignment

    Returns:
        :class:`list` of :class:`~pysam.AlignedSegment`: list of aligned segments

    Note:
        using a higher min_match may improve performance as low quality alignments are rejected more quickly. However
        this may also result in no match being returned when there is no high quality match to be found.
    """
    ref = str(ref)
    if len(ref) < 1 or len(seq) < 1:
        raise AttributeError('cannot overlap on an empty sequence: len(ref)={}, len(seq)={}'.format(len(ref), len(seq)))
    if min_match < 0 or min_match > 1:
        raise AttributeError('min_match must be between 0 and 1')

    if min_overlap_percent <= 0 or min_overlap_percent > 1:
        raise AttributeError('percent must be greater than 0 and up to 1', min_overlap_percent)

    min_overlap = int(round(min_overlap_percent * len(seq), 0))
    # store to improve speed and space (don't need to store all alignments)
    best_score = (0, 0)
    results = []

    putative_start_positions = range(min_overlap - len(seq), len(ref) + len(seq) - min_overlap)
    if min_consecutive_match > 1:
        putative_start_positions = set()
        kmers_checked = {}
        for i in range(0, len(seq) - min_consecutive_match):
            current_kmer = seq[i:i + min_consecutive_match]
            if current_kmer in kmers_checked:
                putative_start_positions.update([p - i for p in kmers_checked[current_kmer]])
                continue
            rp = [m.start() for m in re.finditer(current_kmer, ref)]
            kmers_checked[current_kmer] = rp
            putative_start_positions.update([p - i for p in rp])
    for ref_start in putative_start_positions:
        score = 0
        cigar = []
        mismatches = 0
        length = len(seq)
        for i in range(0, len(seq)):
            if length == 0:
                break
            r = ref_start + i
            if r < 0 or r >= len(ref):  # outside the length of the reference seq
                cigar.append((CIGAR.S, 1))
                length -= 1
                continue
            if DNA_ALPHABET.match(ref[r], seq[i]):
                cigar.append((CIGAR.EQ, 1))
            else:
                cigar.append((CIGAR.X, 1))
                mismatches += 1
                if mismatches / length > 1 - min_match:
                    break
        if length == 0 or mismatches / length > 1 - min_match:
            continue

        cigar = _cigar.join(cigar)
        # end mismatches we set as soft-clipped
        if cigar[0][0] == CIGAR.X:
            cigar[0] = (CIGAR.S, cigar[0][1])
        if cigar[-1][0] == CIGAR.X:
            cigar[-1] = (CIGAR.S, cigar[-1][1])

        qstart = 0 if cigar[0][0] != CIGAR.S else cigar[0][1]

        a = SamRead(
            query_sequence=str(seq),
            reference_start=ref_start + qstart,
            cigar=cigar
        )
        qlen = a.reference_end - a.reference_start
        score = (scoring_function(a), qlen)  # this way for equal identity matches we take the longer alignment
        if qlen < min_overlap:
            continue
        if score >= best_score:
            best_score = score
            results.append((a, score))

    filtered = [x for x, y in results if y == best_score]
    return filtered


def sequenced_strand(read, strand_determining_read=2):
    """
    determines the strand that was sequenced

    Args:
        read (:class:`~pysam.AlignedSegment`): the read being used to determine the strand
        strand_determining_read (int): which read in the read pair is the same as the sequenced strand

    Returns:
        STRAND: the strand that was sequenced

    Raises:
        ValueError: if strand_determining_read is not 1 or 2

    Warning:
        if the input pair is unstranded the information will not be representative of the
        strand sequenced since the assumed convention is not followed
    """
    if read.is_unmapped or not read.is_paired:
        raise ValueError('cannot determine strand if the read is unmapped or unpaired')
    strand = None
    if strand_determining_read == 1:
        if read.is_read1:
            strand = STRAND.NEG if read.is_reverse else STRAND.POS
        else:
            strand = STRAND.NEG if not read.is_reverse else STRAND.POS
    elif strand_determining_read == 2:
        if read.is_read2:
            strand = STRAND.NEG if read.is_reverse else STRAND.POS
        else:
            strand = STRAND.NEG if not read.is_reverse else STRAND.POS
    else:
        raise ValueError('unexpected value. Expected 1 or 2, found:', strand_determining_read)
    return strand


def read_pair_type(read):
    # check if the read pair is in the expected orientation
    """
    assumptions based on illumina pairs: only 4 possible combinations

    Args:
        read (:class:`~pysam.AlignedSegment`): the input read

    Returns:
        READ_PAIR_TYPE: the type of input read pair

    Raises:
        NotImplementedError: for any read that does not fall into the four expected configurations (see below)

    ::

        ++++> <---- is LR same-strand
        ++++> ++++> is LL opposite
        <---- <---- is RR opposite
        <---- ++++> is RL same-strand
    """
    reverse = False
    if read.reference_id > read.next_reference_id or \
            (read.reference_id == read.next_reference_id and read.reference_start > read.next_reference_start):
        reverse = True

    if not read.is_reverse and read.mate_is_reverse:  # LR
        return READ_PAIR_TYPE.RL if reverse else READ_PAIR_TYPE.LR
    elif not read.is_reverse and not read.mate_is_reverse:  # LL opp
        return READ_PAIR_TYPE.LL
    elif read.is_reverse and read.mate_is_reverse:  # RR opp
        return READ_PAIR_TYPE.RR
    elif read.is_reverse and not read.mate_is_reverse:  # RL
        return READ_PAIR_TYPE.LR if reverse else READ_PAIR_TYPE.RL
    else:
        raise NotImplementedError('unexpected orientation for pair')


def orientation_supports_type(read, event_type):
    """
    checks if the orientation is compatible with the type of event

    Args:
        read (:class:`~pysam.AlignedSegment`): a read from the pair
        event_type (SVTYPE): the type of event to check

    Returns:
        bool:
            - ``True`` - the read pair is in the correct orientation for this event type
            - ``False`` - the read is not in the correct orientation
    """
    if event_type == SVTYPE.DEL or event_type == SVTYPE.INS:
        if read_pair_type(read) != READ_PAIR_TYPE.LR:
            return False
    elif event_type == SVTYPE.TRANS:
        if read_pair_type(read) != READ_PAIR_TYPE.LR and \
                read_pair_type(read) != READ_PAIR_TYPE.RL:
            return False
    elif event_type == SVTYPE.ITRANS or event_type == SVTYPE.INV:
        if read_pair_type(read) != READ_PAIR_TYPE.LL and \
                read_pair_type(read) != READ_PAIR_TYPE.RR:
            return False
    elif event_type == SVTYPE.DUP:
        if read_pair_type(read) != READ_PAIR_TYPE.RL:
            return False
    else:
        raise ValueError('unexpected event type', event_type)
    return True


def convert_events_to_softclipping(read, orientation, max_event_size, min_anchor_size=None):
    """
    given an alignment, simplifies the alignment by grouping everything past the first anchor and including the
    first event considered too large and unaligning them turning them into softclipping

    """
    if min_anchor_size is None:
        min_anchor_size = max_event_size

    if orientation == ORIENT.LEFT:
        event_size = 0
        adjusted_cigar = []
        anchor = 0
        for state, count in read.cigar:
            if state == CIGAR.M:
                raise NotImplementedError('match v mismatch must be specified')
            elif anchor < min_anchor_size:
                if state == CIGAR.EQ:
                    anchor += count
            elif state in EVENT_STATES:
                event_size += count
                if event_size > max_event_size:
                    break
            else:
                event_size = 0
            adjusted_cigar.append((state, count))
        if event_size > max_event_size:
            while adjusted_cigar[-1][0] in EVENT_STATES:
                del adjusted_cigar[-1]
            aligned = sum([y for x, y in adjusted_cigar if x in QUERY_ALIGNED_STATES] + [0])
            sc = len(read.query_sequence) - aligned
            adjusted_cigar.append((CIGAR.S, sc))
            read = copy(read)
            read.cigar = adjusted_cigar
    elif orientation == ORIENT.RIGHT:
        # more complicated than left b/c need to also adjust the start position
        event_size = 0
        anchor = 0
        adjusted_cigar = []
        for state, count in read.cigar[::-1]:  # first event from the right
            if state == CIGAR.M:
                raise NotImplementedError('match v mismatch must be specified')
            elif anchor < min_anchor_size:
                if state == CIGAR.EQ:
                    anchor += count
            elif state in EVENT_STATES:
                event_size += count
                if event_size > max_event_size:
                    break
            else:
                event_size = 0
            adjusted_cigar.append((state, count))
        if event_size > max_event_size:
            while adjusted_cigar[-1][0] in EVENT_STATES:
                del adjusted_cigar[-1]
            originally_refaligned = sum([y for x, y in read.cigar if x in REFERENCE_ALIGNED_STATES] + [0])
            refaligned = sum([y for x, y in adjusted_cigar if x in REFERENCE_ALIGNED_STATES] + [0])
            aligned = sum([y for x, y in adjusted_cigar if x in QUERY_ALIGNED_STATES] + [0])
            sc = len(read.query_sequence) - aligned
            adjusted_cigar = [(CIGAR.S, sc)] + adjusted_cigar[::-1]
            read = copy(read)
            read.cigar = adjusted_cigar
            read.reference_start += originally_refaligned - refaligned
    else:
        raise ValueError('orientation must be specified', orientation)
    return read


def consensus_reads(reads, min_consensus_score=0.5, min_consensus_count=3, min_length=0):
    """
    compute a consensus read for all input reads
    """
    ref_intersect = Interval.intersection(*[Interval(r.reference_start, r.reference_end) for r in reads])
    if not ref_intersect:
        raise ValueError('Cannot compute a consensus for reads without a common intersection')
    # left most position of the reference intersection is where all reads should start alignment
    phred_hist = {}  # query_pos => (ref_pos, base, score, cigar_state)
    for read in reads:
        offset = read.reference_start - ref_intersect.start - read.query_alignment_start
        cigar_values = []
        ref_pos = read.reference_start
        for state, freq in read.cigar:
            for _ in range(0, freq):
                if state in _cigar.REFERENCE_ALIGNED_STATES:
                    cigar_values.append((ref_pos, state))
                    ref_pos += 1
                elif state != CIGAR.H:
                    cigar_values.append((None, state))

        for index in range(0, len(read.query_sequence)):
            phred_hist.setdefault(offset + index, []).append((
                read.query_sequence[index],
                read.query_qualities[index] if read.query_qualities else None,
                cigar_values[index][0],
                cigar_values[index][1],
                read  # for tracking the consensus back to each read
            ))
    # finds the consensus (where possible) at all positions
    consensus = []
    for query_pos, values in phred_hist.items():
        if len(values) < min_consensus_count:
            continue
        cons_base, cons_phred, cons_score = consensus_base([(v[0], v[1]) for v in values])
        if cons_score < min_consensus_score:
            continue
        other = {(v[2], v[3]) for v in values if v[0] == cons_base}
        if (None, CIGAR.S) in other and len(other) == 2:
            other = {(None, CIGAR.S)}
        if len(other) != 1:
            continue
        ref_pos, cigar_state = other.pop()
        consensus.append((query_pos, cons_base, cons_phred, ref_pos, cigar_state, [v[4] for v in values if v[0] == cons_base]))

    consecutive_cons = []
    for value in sorted(consensus):
        if not consecutive_cons or consecutive_cons[-1][-1][0] != value[0] - 1:
            consecutive_cons.append([value])
        else:
            consecutive_cons[-1].append(value)
    print('consecutive_cons', [len(v) for v in consecutive_cons])

    read_support = {}  # read => {cons read => score}
    for curr in sorted(consecutive_cons, key=lambda x: len(x), reverse=True):
        if len(curr) < min_length:
            continue
        if not {c[4] for c in curr} - {CIGAR.S, CIGAR.H}:  # nothing aligned
            continue
        cigar = []
        for _, _, _, ref_pos, cigar_state, _ in curr:
            if not cigar or cigar[-1][1] != cigar_state:
                cigar.append((ref_pos, cigar_state, 1))
            else:
                cigar[-1] = (cigar[-1][0], cigar_state, cigar[-1][2] + 1)
        print(cigar)
        for i, (_, state, _) in enumerate(cigar):
            if state == CIGAR.S and i > 0 and i < len(cigar) - 1:
                prev_align_count = 0
                for _, prev_state, prev_freq in cigar[:i:-1]:
                    if prev_state not in {CIGAR.H, CIGAR.S}:
                        prev_align_count += prev_freq
                    else:
                        break
                next_align_count = 0
                for _, next_state, next_freq in cigar[i + 1:]:
                    if prev_state not in {CIGAR.H, CIGAR.S}:
                        next_align_count += next_freq
                    else:
                        break
                if prev_align_count >= next_align_count:
                    for i in range(0, i):
                        cigar[i] = (cigar[i][0], CIGAR.S, cigar[i][2])
                else:
                    for i in range(i + 1, len(cigar)):
                        cigar[i] = (cigar[i][0], CIGAR.S, cigar[i][2])
        print(cigar)
        cons_read = SamRead(
            reference_start=min([c[0] for c in cigar if c[1] in _cigar.REFERENCE_ALIGNED_STATES]),
            query_sequence=''.join([c[1] for c in curr]),
            cigar=_cigar.join([(c[1], c[2]) for c in cigar]),
            query_qualities=[c[2] for c in curr]
        )
        print('cons_read', _cigar.join([(c[1], c[2]) for c in cigar]))
        print(cons_read)
        for reads in [v[5] for v in curr]:
            for read in reads:
                read_support.setdefault(read, {})
                read_support[read][cons_read] = read_support[read].get(cons_read, 0) + 1
    print()
    cons_to_reads = {}
    # now assign each read to its best scoring consensus
    for read in read_support:
        best_read, best_score = max(read_support[read].items(), key=lambda x: x[1])
        best_score = best_score / min(len(best_read.query_sequence), len(read.query_sequence))
        if best_score < min_consensus_score:
            print('\nscore too low')
            print(read)
            print(best_read)
            print(best_score)
            continue
        cons_to_reads.setdefault(best_read, []).append(read)
    for read, support in cons_to_reads.items():
        print(read)
        for r in support:
            print('supp by', r)


def consensus_base(bases):
    """
    For a list of tuples representing bases and phred scores computes the consensus base, weighted by phred score

    Returns:
        tuple:
            - str: the consensus base
            - int: the new phred (average phred)
            - float: the score for the consensus (0-1)
    """
    total_phred = 0
    probs = {}
    for base, phred_qual in bases:
        probs.setdefault(base.upper(), []).append(phred_qual)
        total_phred += phred_qual
    max_base, max_phred = min([(k, v) for k, v in probs.items()], key=lambda x: (-1 * sum(x[1]), x[0]))
    return max_base, int(round(sum(max_phred) / len(max_phred), 0)), sum(max_phred) / total_phred
