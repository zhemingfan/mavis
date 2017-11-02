import json
import os
import re
import time
import warnings

from .constants import DEFAULTS
from .genomic import UsTranscript
from .variant import annotate_events, choose_more_annotated, choose_transcripts_by_priority, determine_prime, flatten_fusion_transcript, flatten_fusion_translation
from ..constants import COLUMNS, PRIME, PROTOCOL, sort_columns
from ..error import DrawingFitError, NotSpecifiedError
from ..illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from ..illustrate.constants import DiagramSettings
from ..illustrate.diagram import draw_sv_summary_diagram
from ..util import generate_complete_stamp, log, mkdirp, read_inputs


ACCEPTED_FILTERS = {
    'choose_more_annotated': choose_more_annotated,
    'choose_transcripts_by_priority': choose_transcripts_by_priority
}


def draw(drawing_config, ann, reference_genome, template_metadata, drawings_directory):
    """
    produces the svg diagram and json legend for a given annotation
    """
    drawing = None
    legend = None
    retry_count = 0
    draw_fusion_transcript = True
    draw_reference_transcripts = True
    initial_width = drawing_config.width
    while drawing is None:  # continue if drawing error and increase width
        try:
            canvas, legend_json = draw_sv_summary_diagram(
                drawing_config, ann, reference_genome=reference_genome,
                templates=template_metadata,
                draw_fusion_transcript=draw_fusion_transcript,
                draw_reference_transcripts=draw_reference_transcripts
            )

            gene_aliases1 = 'NA'
            gene_aliases2 = 'NA'
            try:
                if ann.transcript1.gene.aliases:
                    gene_aliases1 = '-'.join(ann.transcript1.gene.aliases)
                if ann.transcript1.is_best_transcript:
                    gene_aliases1 = 'b-' + gene_aliases1
            except AttributeError:
                pass
            try:
                if ann.transcript2.gene.aliases:
                    gene_aliases2 = '-'.join(ann.transcript2.gene.aliases)
                if ann.transcript2.is_best_transcript:
                    gene_aliases2 = 'b-' + gene_aliases2
            except AttributeError:
                pass
            try:
                if determine_prime(ann.transcript1, ann.break1) == PRIME.THREE:
                    gene_aliases1, gene_aliases2 = gene_aliases2, gene_aliases1
            except NotSpecifiedError:
                pass

            name = 'mavis_{}-chr{}_chr{}-{}_{}'.format(
                ann.annotation_id, ann.break1.chr, ann.break2.chr, gene_aliases1, gene_aliases2
            )

            drawing = os.path.join(drawings_directory, name + '.svg')
            legend = os.path.join(drawings_directory, name + '.legend.json')
            log('generating svg:', drawing, time_stamp=False)
            canvas.saveas(drawing)

            log('generating legend:', legend, time_stamp=False)
            with open(legend, 'w') as fh:
                json.dump(legend_json, fh)
            break
        except DrawingFitError as err:
            drawing_config.width += drawing_config.drawing_width_iter_increase
            log('extending width by', drawing_config.drawing_width_iter_increase, 'to', drawing_config.width, time_stamp=False)
            retry_count += 1
            if retry_count > drawing_config.max_drawing_retries:
                if draw_fusion_transcript and draw_reference_transcripts:
                    log('restricting to gene-level only', time_stamp=False)
                    draw_fusion_transcript = False
                    draw_reference_transcripts = False
                    drawing_config.width = initial_width
                    retry_count = 0
                else:
                    warnings.warn(str(err))
                    drawing = True
    drawing_config.width = initial_width  # reset the width
    return drawing, legend


def main(
    inputs, output, library, protocol,
    reference_genome, annotations, template_metadata,
    min_domain_mapping_match=DEFAULTS.min_domain_mapping_match,
    min_orf_size=DEFAULTS.min_orf_size,
    max_orf_cap=DEFAULTS.max_orf_cap,
    annotation_filters=DEFAULTS.annotation_filters,
    start_time=int(time.time()),
    draw_fusions_only=DEFAULTS.draw_fusions_only,
    **kwargs
):
    """
    Args:
        inputs (:class:`List` of :class:`str`): list of input files to read
        output (str): path to the output directory
        reference_genome (object): see :func:`~mavis.annotate.file_io.load_reference_genome`
        annotations(object): see :func:`~mavis.annotate.file_io.load_reference_genes`
        template_metadata (object): see :func:`~mavis.annotate.file_io.load_templates`
        min_domain_mapping_match (float): min mapping match percent (0-1) to count a domain as mapped
        min_orf_size (int): minimum size of an :term:`open reading frame` to keep as a putative translation
        max_orf_cap (int): the maximum number of :term:`open reading frame` s to collect for any given event
    """
    drawings_directory = os.path.join(output, 'drawings')
    tabbed_output_file = os.path.join(output, 'annotations.tab')
    fa_output_file = os.path.join(output, 'annotations.fusion-cdna.fa')

    annotation_filters = [] if not annotation_filters else annotation_filters.split(',')
    annotation_filters = [ACCEPTED_FILTERS[a] for a in annotation_filters]

    mkdirp(drawings_directory)
    # test that the sequence makes sense for a random transcript
    bpps = read_inputs(
        inputs, in_={COLUMNS.protocol: PROTOCOL.values()},
        add_default={
            COLUMNS.protocol: protocol,
            COLUMNS.library: library,
            COLUMNS.stranded: False
        },
        require=[COLUMNS.protocol, COLUMNS.library],
        expand_strand=False, expand_orient=True, expand_svtype=True
    )
    log('read {} breakpoint pairs'.format(len(bpps)))

    annotations = annotate_events(
        bpps,
        reference_genome=reference_genome, annotations=annotations,
        min_orf_size=min_orf_size, min_domain_mapping_match=min_domain_mapping_match,
        max_orf_cap=max_orf_cap,
        log=log,
        filters=annotation_filters
    )

    fa_sequence_names = set()
    # now try generating the svg
    drawing_config = DiagramSettings(**{k: v for k, v in kwargs.items() if k in ILLUSTRATION_DEFAULTS})

    header_req = {
        COLUMNS.break1_strand,
        COLUMNS.break2_strand,
        COLUMNS.fusion_sequence_fasta_file,
        COLUMNS.fusion_splicing_pattern,
        COLUMNS.fusion_cdna_coding_start,
        COLUMNS.fusion_cdna_coding_end,
        COLUMNS.fusion_sequence_fasta_id,
        COLUMNS.fusion_mapped_domains,
        COLUMNS.exon_first_3prime,
        COLUMNS.exon_last_5prime,
        COLUMNS.annotation_id,
        COLUMNS.annotation_figure,
        COLUMNS.annotation_figure_legend,
        COLUMNS.cdna_synon,
        COLUMNS.protein_synon
    }
    header = None
    log('opening for write:', tabbed_output_file)
    tabbed_fh = open(tabbed_output_file, 'w')
    log('opening for write:', fa_output_file)
    fasta_fh = open(fa_output_file, 'w')

    try:
        total = len(annotations)
        for i, ann in enumerate(annotations):
            ann_row = ann.flatten()
            ann_row[COLUMNS.fusion_sequence_fasta_file] = fa_output_file
            if header is None:
                header_req.update(ann_row.keys())
                header = sort_columns(header_req)
                tabbed_fh.write('\t'.join([str(c) for c in header]) + '\n')
            log(
                '({} of {}) current annotation'.format(i + 1, total),
                ann.annotation_id, ann.transcript1, ann.transcript2, ann.event_type)
            # get the reference sequences for either transcript
            ref_cdna_seq = {}
            ref_protein_seq = {}

            for ust in [x for x in [ann.transcript1, ann.transcript2] if isinstance(x, UsTranscript)]:
                name = ust.name
                for spl_tx in ust.spliced_transcripts:
                    ref_cdna_seq.setdefault(spl_tx.get_seq(reference_genome), set()).add(name)
                    for translation in spl_tx.translations:
                        ref_protein_seq.setdefault(translation.get_aa_seq(reference_genome), set()).add(name)

            # try building the fusion product
            rows = []
            # add fusion information to the current ann_row
            for spl_fusion_tx in [] if not ann.fusion else ann.fusion.transcripts:
                fusion_fa_id = '{}_{}'.format(ann.annotation_id, spl_fusion_tx.splicing_pattern.splice_type)
                fusion_fa_id = re.sub(r'\s', '-', fusion_fa_id)
                if fusion_fa_id in fa_sequence_names:
                    raise AssertionError('should not be duplicate fa sequence ids', fusion_fa_id)
                seq = ann.fusion.get_cdna_seq(spl_fusion_tx.splicing_pattern)
                fasta_fh.write('> {}\n{}\n'.format(fusion_fa_id, seq))
                cdna_synon = ';'.join(sorted(list(ref_cdna_seq.get(seq, set()))))

                temp_row = {}
                temp_row.update(ann_row)
                temp_row.update(flatten_fusion_transcript(spl_fusion_tx))
                temp_row[COLUMNS.fusion_sequence_fasta_id] = fusion_fa_id
                temp_row[COLUMNS.cdna_synon] = cdna_synon if cdna_synon else None
                if spl_fusion_tx.translations:
                    # duplicate the ann_row for each translation
                    for fusion_translation in spl_fusion_tx.translations:
                        nrow = dict()
                        nrow.update(ann_row)
                        nrow.update(temp_row)
                        aa_seq = fusion_translation.get_aa_seq()
                        protein_synon = ';'.join(sorted(list(ref_protein_seq.get(aa_seq, set()))))
                        nrow[COLUMNS.protein_synon] = protein_synon if protein_synon else None
                        # select the exon
                        nrow.update(flatten_fusion_translation(fusion_translation))
                        rows.append(nrow)
                else:
                    temp_row.update(ann_row)
                    rows.append(temp_row)
            # draw the annotation and add the path to all applicable rows (one drawing for multiple annotations)
            if ann.fusion or not draw_fusions_only:
                drawing, legend = draw(drawing_config, ann, reference_genome, template_metadata, drawings_directory)
                for row in rows + [ann_row]:
                    row[COLUMNS.annotation_figure] = drawing
                    row[COLUMNS.annotation_figure_legend] = legend
            if not rows:
                rows = [ann_row]
            for row in rows:
                tabbed_fh.write('\t'.join([str(row.get(k, None)) for k in header]) + '\n')
        generate_complete_stamp(output, log, start_time=start_time)
    finally:
        log('closing:', tabbed_output_file)
        tabbed_fh.close()
        log('closing:', fa_output_file)
        fasta_fh.close()
