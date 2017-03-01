from .constants import *
from .error import *


class Interval:
    """
    """
    def __init__(self, start, end=None, freq=1, number_type=None):
        """
        Args:
            start (int): the start of the interval (inclusive)
            end (int): the end of the interval (inclusive)
            freq (int): the frequency or weight of the interval
        """
        self.start = start
        self.end = end if end is not None else start
        
        if number_type is None:
            if int(self.start) != float(self.start) or int(self.end) != float(self.end) \
                    or type(self.start) == float or type(self.end) == float:
                number_type = float
            else:
                number_type = int
        self.number_type = number_type

        self.start = self.number_type(self.start)
        self.end = self.number_type(self.end)
        if self.start > self.end:
            raise AttributeError('interval start > end is not allowed', self.start, self.end)
        self.freq = int(freq)
        if self.freq <= 0:
            raise AttributeError('Interval frequency must be a natural number')

    def __sub__(self, other):  # difference
        """the difference of two intervals

        Example:
            >>> Interval(1, 10) - Interval(5, 50)
            [Interval(1, 4)]
            >>> Interval(1, 2) - Interval(10, 11)
            [Interval(1, 2)]
            >>> Interval(1, 2) - Interval(-1, 10)
            []
        """
        if Interval.overlaps(self, other):
            if other[0] <= self[0]:
                if other[1] >= self[1]:
                    return []
                else:
                    return [Interval(other[1] + 1, self[1])]
            elif other[1] >= self[1]:
                return [Interval(self[0], other[0] - 1)]
            else:
                return [Interval(self[0], other[0] - 1), Interval(other[1] + 1, self[1])]
        else:
            return [Interval(self[0], self[1])]


    def __and__(self, other):  # intersection
        """the intersection of two intervals

        Example:
            >>> Interval(1, 10) & Interval(5, 50)
            Interval(5, 10)
            >>> Interval(1, 2) & Interval(10, 11)
            None
        """
        return Interval.intersection(self, other)

    def __or__(self, other):  # union
        """the union of two intervals

        Example:
            >>> Interval(1, 10) | Interval(5, 50)
            Interval(1, 50)
            >>> Interval(1, 2) | Interval(10, 11)
            Interval(1, 11)
        """
        return Interval.union(self, other)

    def __xor__(self, other):
        return (self - other) + (other - self)

    def __getitem__(self, index):
        try:
            index = int(index)
        except ValueError:
            raise IndexError('index input accessor must be an integer', index)
        if index == 0:
            return self.start
        elif index == 1:
            return self.end
        raise IndexError(
            'index input accessor is out of bounds: 1 or 2 only', index)

    @classmethod
    def overlaps(cls, first, other):
        """
        checks if two intervals have any portion of their given ranges in common

        Args:
            first (Interval): an interval to be compared
            other (Interval): an interval to be compared

        Example:
            >>> Interval.overlaps(Interval(1, 4), Interval(5, 7))
            False
            >>> Interval.overlaps(Interval(1, 10), Interval(10, 11))
            True
            >>> Interval.overlaps((1, 10), (10, 11))
            True
        """
        if first[1] < other[0]:
            return False
        elif first[0] > other[1]:
            return False
        else:
            return True

    def __len__(self):
        """
        the length of the interval

        Example:
            >>> len(Interval(1, 11))
            12

        Warning:
            only works for integer intervals
        """
        return Interval.length(self)

    def length(self):
        try:
            if self.number_type == float:
                return self[1] - self[0]
        except AttributeError:
            pass
        return self[1] - self[0] + 1

    def __lt__(self, other):
        if self[0] < other[0]:
            return True
        elif self[0] == other[0] and self[1] < other[1]:
            return True
        return False

    def __gt__(self, other):
        if self[1] > other[1]:
            return True
        elif self[1] == other[1] and self[0] > other[0]:
            return True
        return False

    def __repr__(self):
        cls = self.__class__.__name__
        number_type = ''
        if self.number_type != int:
            number_type = ', type={}'.format(self.number_type)
        if self.freq != 1:
            return '{}({}, {}, freq={}{})'.format(cls, self.start, self.end, self.freq, number_type)
        else:
            return '{}({}, {}{})'.format(cls, self.start, self.end, number_type)

    @property
    def center(self):
        """
        the middle of the interval

        Example:
            >>> Interval(1, 10).center
            5.5
            >>> Interval(1, 11).center
            6
        """
        return (self[1] + self[0]) / 2

    def __eq__(self, other):
        if self[0] != other[0] or self[1] != other[1]:
            return False
        return True

    def __contains__(self, other):
        if other[0] >= self[0] and other[1] <= self[1]:
            return True
        return False

    @classmethod
    def dist(cls, first, other):
        """returns the minimum distance between intervals

        Example:
            >>> Interval.dist((1, 4), (5, 7))
            -1
            >>> Interval.dist((5, 7), (1, 4))
            1
            >>> Interval.dist((5, 8), (7, 9))
            0
        """
        if first[1] < other[0]:
            return first[1] - other[0]
        elif first[0] > other[1]:
            return first[0] - other[1]
        else:
            return 0

    def __hash__(self):
        return hash((self[0], self[1], self.freq))

    # @classmethod
    # def weighted_mean_ci(cls, *intervals):
    #     """
    #     Calculates the weighted mean of a set of intervals. The weighting is inversely proportional to
    #     the size of the input interval. The length of the final interval is the weight mean length of 
    #     the input intervals also weighted by length

    #     Args:
    #         intervals (Interval): a list of intervals

    #     Returns:
    #         Interval: the weighted mean interval of the input intervals

    #     Raises:
    #         AttributeError: if the input list is empty

    #     Example:
    #         >>> Interval.weighted_mean((1, 2), (1, 9), (2, 10))
    #         Interval(1, 4)
    #         >>> Interval.weighted_mean((1, 1), (10, 10))
    #         Interval(6)
    #     """
    #     centers = []
    #     weights = []
    #     lengths = []
    #     number_type = int
    #     if len(intervals) == 0:
    #         raise AttributeError('cannot compute the weighted mean interval of an empty set of intervals')
    #     for i in intervals:
    #         try:
    #             if i.number_type == float:
    #                 number_type = float
    #         except AttributeError:
    #             pass
    #         if not isinstance(i, Interval):
    #             i = Interval(i[0], i[1])
    #         for temp in range(0, i.freq):
    #             centers.append(i.center)
    #             weights.append(1 / i.length())
    #             lengths.append(i.length())

    #     center = np.average(centers, weights=weights)
    #     size = max([np.average(lengths) - 1, 1]) if number_type == int else np.average(lengths)
    #     return Interval(round(center - size / 2, 0), round(center + size / 2, 0))

    @classmethod
    def position_in_range(cls, segments, pos):
        if len(segments) == 0:
            raise AttributeError('cannot compute on an empty list')

        num = 0
        found_inbetween_segment = False

        segments = sorted(segments, key=lambda x: (x[0], x[1]))

        while num < len(segments):
            current = segments[num]

            if pos[0] >= current[0] \
                    and pos[1] <= current[1]:
                # pos range is fully contained in the current segment
                break
            elif num == 0:  # first segment
                if pos[1] < current[0]:
                    # before the first segment
                    found_inbetween_segment = True
                    break
            else:
                # check the previous segment
                previous = segments[num - 1]
                if pos[0] > previous[1] and pos[1] < current[0]:
                    found_inbetween_segment = True
                    break
            num += 1
        return num, found_inbetween_segment

    @classmethod
    def convert_pos(cls, mapping, pos, forward_to_reverse=None):
        i = cls.convert_ratioed_pos(mapping, pos, forward_to_reverse)
        if i.forward_to_reverse:
            return i.end
        else:
            return i.start

    @classmethod
    def convert_ratioed_pos(cls, mapping, pos, forward_to_reverse=None):
        """ convert any given position given a mapping of intervals to another range

        Args:
            mapping (:class:`dict` of :class:`Interval` and :class:`Interval`): a mapping of a set of continuous intervals
            pos (int): a position in the first coordinate system

        Returns:
            int: the position in the alternate coordinate system given the input mapping

        Raises:
            AttributeError: if the input position is outside the set of input segments
            IndexError: if the input position cannot be converted to the output system

        Example:
            >>> mapping = {(1, 10): (101, 110), (11, 20): (555, 564)}
            >>> Interval.convert_pos(mapping, 5)
            5
            >>> Interval.convert_pos(mapping, 15)
            559
        """
        if len(mapping.keys()) < 2 and forward_to_reverse is None:
            raise AttributeError(
                'mapping is insufficient to determine orientation')

        # order the input intervals
        input_intervals = sorted(mapping.keys())
        mapped_to_intervals = [mapping[i] for i in input_intervals]

        # input check interval ranges are non-overlapping and increasing/decreasing
        for i in range(0, len(input_intervals)):
            if i > 0:
                if Interval.overlaps(input_intervals[i - 1], input_intervals[i]):
                    raise AttributeError(
                        'input intervals cannot be overlapping',
                        input_intervals[i], input_intervals[i - 1]
                    )
                """if Interval.overlaps(mapped_to_intervals[i - 1], mapped_to_intervals[i]):
                    raise AttributeError(
                        'mapped_to intervals cannot be overlapping',
                        mapped_to_intervals[i], mapped_to_intervals[i - 1]
                    )"""
                if mapped_to_intervals[i][0] > mapped_to_intervals[i - 1][1]:
                    if forward_to_reverse is None:
                        forward_to_reverse = False
                    elif forward_to_reverse:
                        raise AttributeError('direction of mapped intervals is not consistent')
                elif mapped_to_intervals[i][1] < mapped_to_intervals[i - 1][0]:
                    if forward_to_reverse is None:
                        forward_to_reverse = True
                    elif not forward_to_reverse:
                        raise AttributeError('direction of mapped intervals is not consistent')

        i, previous_flag = Interval.position_in_range(
            input_intervals, (pos, pos))  # get the input position
        if i == len(input_intervals):
            curr = input_intervals[i - 1]
            if not forward_to_reverse:
                raise IndexError(pos, 'is outside mapped range', mapping)
            else:
                raise IndexError(pos, 'is outside mapped range', mapping)
        elif previous_flag:
            curr = input_intervals[i]
            if i == 0:
                if not forward_to_reverse:
                    raise IndexError(pos, 'is outside mapped range', mapping)
                else:
                    raise IndexError(pos, 'is outside mapped range', mapping)
            else:  # between two segments
                prev = input_intervals[i - 1]
                if not forward_to_reverse:
                    raise IndexError(pos, 'is outside mapped range', mapping)
                else:
                    raise IndexError(pos, 'is outside mapped range', mapping)
        else:
            # fell into a mapped region
            curr = input_intervals[i]
            nexxt = mapping[curr]
            if curr[1] - curr[0] == 0:
                i = Interval(nexxt[0], nexxt[1])
            else:
                ratio = (nexxt[1] - nexxt[0]) / (curr[1] - curr[0])
                shift = round((pos - curr[0]) * ratio, 0)
                shift2 = round((pos - curr[0]) * ratio + ratio, 0)
                number_type = int if ratio == 1 else float
                if forward_to_reverse:
                    i = Interval(nexxt[1] - shift2, nexxt[1] - shift, number_type=number_type)
                else:
                    i = Interval(nexxt[0] + shift, nexxt[0] + shift2, number_type=number_type)
            setattr(i, 'forward_to_reverse', forward_to_reverse)
            return i

    @classmethod
    def union(cls, *intervals):
        """
        returns the union of the set of input intervals

        Example:
            >>> Interval.union((1, 2), (4, 6), (4, 9), (20, 21))
            Interval(1, 21)
        """
        if len(intervals) < 1:
            raise AttributeError('cannot compute the union of an empty set of intervals')
        return Interval(min([i[0] for i in intervals]), max([i[1] for i in intervals]))

    @classmethod
    def intersection(cls, *intervals):
        """
        returns None if there is no intersection

        Example:
            >>> Interval.intersection((1, 10), (2, 8), (7, 15))
            Interval(7, 8)
            >>> Interval.intersection((1, 2), (5, 9))
            None
        """
        if len(intervals) < 1:
            raise AttributeError('cannot compute the intersection of an empty set of intervals')
        low = max([i[0] for i in intervals])
        high = min([i[1] for i in intervals])
        if low > high:
            return None
        return Interval(low, high)

    @classmethod
    def min_nonoverlapping(cls, *intervals):
        """
        for a list of intervals, orders them and merges any overlap to return a list of non-overlapping intervals
        O(nlogn)

        Example:
            >>> Interval.min_nonoverlapping((1, 10), (7, 8), (6, 14), (17, 20))
            [Interval(1, 14), Interval(17, 20)]
        """
        if len(intervals) == 0:
            return []
        intervals = sorted(list(intervals), key=lambda x: (x[0], x[1]))
        new_intervals = [Interval(intervals[0][0], intervals[0][1])]
        for i in intervals[1:]:
            if Interval.overlaps(new_intervals[-1], i):
                new_intervals[-1] = new_intervals[-1] | i
            else:
                new_intervals.append(Interval(i[0], i[1]))
        return new_intervals