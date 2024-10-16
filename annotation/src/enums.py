import enum


class CNVType(enum.StrEnum):
    LOSS = "loss"
    GAIN = "gain"


class OverlapType(enum.Enum):
    all = 0  # return all types of overlap
    span_whole_only = 1  # return only ranges, that spans the full query range
    partial_start = 2  # return only ranges, that overlap start of the range (and not end)
    partial_end = 3  # return only ranges, that overlap end of the range (and not start)
    partial_both = 4  # return only ranges, that overlap end xor start (and not the second)
    inside_only = 5  # return only ranges, that are completely inside the full query range
