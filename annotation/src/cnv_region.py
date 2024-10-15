import re
from dataclasses import dataclass
from functools import cached_property
from annotation.src import enums, constants, core, exceptions


@dataclass
class CNVRegion:
    chr: str
    start: int
    end: int
    cnv_type: enums.CNVType

    def __post_init__(self):
        if self.start > self.end:
            raise exceptions.InvalidRegionRangeError(self.start, self.end)

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def name(self) -> str:
        return f"{self.chr}_{self.start}_{self.end}_{self.cnv_type}"

    @cached_property
    def cytogenetic_position(self) -> str:
        return self._infer_cytogenetic_position()

    @classmethod
    def build_from_string(cls, input_str: str) -> "CNVRegion":
        """
        Build CNVRegion object from string in the format "chr1:10000-20000/del"
        """
        match = re.match(core.REGION_REGEX, input_str)
        if not match:
            raise exceptions.InvalidRegionFormatError(input_str, core.REGION_REGEX)
        if match.group(1) not in constants.ALLOWED_CHROMOSOMES:
            raise exceptions.InvalidRegionChromosomeError(match.group(1), constants.ALLOWED_CHROMOSOMES)

        cnv_type_map = {"del": "loss", "dup": "gain", "gain": "gain", "loss": "loss", "aoh": "loss"}
        if match.group(4).lower() not in cnv_type_map.keys():
            raise exceptions.InvalidRegionCNVTypeError(match.group(4), list(cnv_type_map.keys()))

        return cls(
            chr=match.group(1),
            start=int(match.group(2)),
            end=int(match.group(3)),
            cnv_type=enums.CNVType(cnv_type_map[match.group(4).lower()]),
        )


    def _get_cytobands_overlap(self) -> list[str]:
        import pandas as pd
        cytobands = pd.read_csv(core.CYTOBANDS_FILEPATH, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'color'])
        cyto_overlap = (cytobands['chrom'] == self.chr) & (cytobands['end'] >= self.start) & (cytobands['start'] <= self.end)
        return list(cytobands[cyto_overlap]['name'])


    def _infer_cytogenetic_position(self) -> str:
        """
        Get cytoband description/range as a human-readable string.
        :return: str - cytobands description of the range
        """
        all_cytobands = self._get_cytobands_overlap()

        if len(all_cytobands) == 0:
            raise exceptions.UnknownCytobandError(self.chr, self.start, self.end)

        if len(all_cytobands) == 1:
            return all_cytobands[0]

        return f'{all_cytobands[0]}-{all_cytobands[-1]}'
