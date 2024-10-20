class GenovisioAnnotationError(Exception):
    """Base class for exceptions in this module."""

    pass


class UnknownCytobandError(GenovisioAnnotationError):
    """Exception raised when cytoband could not be determined."""

    def __init__(self, chrom: str, start: int, end: int):
        message = f"Could not determine cytobad for {chrom}:{start}-{end}"
        super().__init__(message)


class InvalidRegionRangeError(GenovisioAnnotationError):
    """Exception raised when region range is invalid."""

    def __init__(self, start: int, end: int):
        message = f"Start position ({start}) must be less than or equal to end position ({end})"
        super().__init__(message)


class InvalidRegionChromosomeError(GenovisioAnnotationError):
    """Exception raised when chromosome name of region is invalid."""

    def __init__(self, parsed_chromosome: str, allowed_chromosomes: list[str]):
        message = f"{parsed_chromosome} is not among {allowed_chromosomes=}"
        super().__init__(message)


class InvalidRegionFormatError(GenovisioAnnotationError):
    """Exception raised when region format is invalid by the regex"""

    def __init__(self, input_string: str, regex: str):
        message = f"{input_string=} does not match the {regex=}, e.g. 'chr1:10000-20000/del'"
        super().__init__(message)


class InvalidRegionCNVTypeError(GenovisioAnnotationError):
    """Exception is raised when CNV type is not in the allowed list."""

    def __init__(self, cnv_type: str, allowed_types: list[str]):
        message = f"{cnv_type=} is not among {allowed_types=}"
        super().__init__(message)


class CNVTypeNormalizationError(GenovisioAnnotationError):
    """Exception is raised when normalizing CNV type fails."""

    def __init__(self, cnv_type: str):
        message = f"{cnv_type=} could not be normalized into a valid CNVType"
        super().__init__(message)


class HighRiskForDuplicationError(GenovisioAnnotationError):
    """Exception is raised when the region is duplication and high risk genes are requested"""

    def __init__(self) -> None:
        message = "Only losses should be considered for high risk genes"
        super().__init__(message)


class UnknownOverlapTypeError(GenovisioAnnotationError):
    """Exception is raised when the overlap type is not in the allowed list."""

    def __init__(self, overlap_type: str, allowed_types: list[str]):
        message = f"{overlap_type=} is not among {allowed_types=}"
        super().__init__(message)


class InvalidJSONAnnotationError(GenovisioAnnotationError):
    """Exception is raised when the JSON is not a valid annotation."""

    def __init__(self, json_file: str, err: str):
        message = f"{json_file=} does not contain a valid annotation of this version. KeyError: {err}"
        super().__init__(message)
