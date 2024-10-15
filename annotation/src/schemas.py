import enum
import gzip
import json
import math
import os
from dataclasses import dataclass
from functools import cached_property
from typing import Any, TypedDict

from annotation.src import core, enums, exceptions


class GenesDBGeneTypes(enum.StrEnum):
    """Types of genes as stored in the GenesDB database."""

    PROTEIN_CODING = "protein_coding"
    PSEUDOGENE = "pseudogene"
    LINC_RNA = "lncRNA"
    R_RNA = "rRNA"
    S_NRNA = "snRNA"
    MIRNA = "miRNA"


class GenesDBGeneTypesCounter(TypedDict):
    protein_coding: int
    pseudogenes: int
    lncrna: int
    rrna: int
    snrna: int
    mirna: int
    other: int


class RegulatoryTypes(enum.StrEnum):
    """Types of regulatory elements as stored in the database."""

    ENHANCER = "enhancer"
    PROMOTER = "promoter"
    OPEN_CHROMATIN_REGION = "open_chromatin_region"
    CTCF_BINDING_SITE = "CTCF_binding_site"
    TF_BINDING_SITE = "TF_binding_site"
    CURATED = "regulatory_curated"
    FLANKING_REGION = "flanking_region"
    SILENCER = "silencer"
    TRANSCRIPTIONAL_CIS_REGULATORY_REGION = "transcriptional_cis_regulatory_region"
    DNASE_I_HYPERSENSITIVE_SITE = "DNase_I_hypersensitive_site"
    ENHANCER_BLOCKING_ELEMENT = "enhancer_blocking_element"
    TATA_BOX = "TATA_box"


class RegulatoryTypesCounter(TypedDict):
    enhancer: int
    promoter: int
    open_chromatin_region: int
    CTCF_binding_site: int
    TF_binding_site: int
    curated: int
    flanking_region: int
    silencer: int
    transcriptional_cis_regulatory_region: int
    DNase_I_hypersensitive_site: int
    enhancer_blocking_element: int
    TATA_box: int
    other: int


class TranscriptRegion(TypedDict):
    length: int
    identifier: str
    cds_overlaps: list[int]
    five_prime_utr_overlaps: list[int]
    three_prime_utr_overlaps: list[int]

    flag_five_inside: bool
    flag_three_inside: bool
    flag_contained: bool


class CommonVariabilityRegion(TypedDict):
    population: str
    frequency: float


def _normalize_type(sv_type: str) -> enums.CNVType:
    if str(sv_type).upper() in ["DUP", "DUPLICATION", "GAIN", "INSERTION"]:
        return enums.CNVType.GAIN
    if str(sv_type).upper() in ["DEL", "DELETION", "LOSS"]:
        return enums.CNVType.LOSS
    raise exceptions.CNVTypeNormalizationError(sv_type)


@dataclass
class CNVRegionAnnotation:
    chr: str
    start: int
    end: int
    cnv_type: enums.CNVType
    length: int
    name: str
    cytogenetic_position: str

    @property
    def is_duplication(self) -> bool:
        return self.cnv_type == enums.CNVType.GAIN

    @property
    def genomic_coord(self) -> str:
        return f"{self.chr}:{self.start}-{self.end}"

    def matches_cnv_type_or_both(self, benign_cnv_type: str) -> bool:
        try:
            return _normalize_type(benign_cnv_type) == self.cnv_type
        except exceptions.CNVTypeNormalizationError:
            if benign_cnv_type.upper() in ["GAIN+LOSS", "LOSS+GAIN"]:
                return True
            else:
                return False  # TODO or raise error?

    def has_start_in_region(self, target_start: int, target_end: int) -> bool:
        """Check if the CNV region has its start in the target region"""
        return target_start <= self.start <= target_end

    def has_end_in_region(self, target_start: int, target_end: int) -> bool:
        """Check if the CNV region has its end in the target region"""
        return target_start <= self.end <= target_end

    def is_overlapping(self, target_start: int, target_end: int, overlap: enums.OverlapType) -> bool:
        """Check if the CNV region overlaps with another region in the specified manner"""
        if overlap == enums.OverlapType.all:
            return target_start < self.end and target_end > self.start
        elif overlap == enums.OverlapType.span_whole_only:
            return target_start <= self.start and target_end >= self.end
        elif overlap == enums.OverlapType.partial_start:
            return target_start <= self.start < target_end <= self.end
        elif overlap == enums.OverlapType.partial_end:
            return self.start < target_start <= self.end <= target_end
        elif overlap == enums.OverlapType.partial_both:
            return (target_start <= self.start < target_end < self.end) or (
                self.start < target_start < self.end <= target_end
            )
        elif overlap == enums.OverlapType.inside_only:
            return self.start <= target_start and target_end <= self.end
        else:
            raise ValueError(f"Unknown overlap type: {overlap}")

    def get_overlap_with_region(self, target_start: int, target_end: int) -> int:
        """Get the number of overlapping bases with another interval."""
        return max(0, min(target_end, self.end) - max(target_start, self.start))

    def is_position_inside(self, position: int) -> bool:
        """Check if the position is inside the CNV region"""
        return self.start <= position <= self.end


@dataclass(frozen=True)
class Annotation:
    cnv: CNVRegionAnnotation
    _benign_cnv: list[dict[str, Any]]
    _benign_cnv_gs_inner: list[dict[str, Any]]
    _benign_cnv_gs_outer: list[dict[str, Any]]
    _regulatory: list[dict[str, Any]]
    _gnomad: list[dict[str, Any]]
    _hi_gene: list[dict[str, Any]]
    _hi_region: list[dict[str, Any]]
    _genes: list[dict[str, Any]]

    @classmethod
    def load_from_json(cls, json_file: str) -> "Annotation":
        if not os.path.exists(json_file):
            raise FileNotFoundError(f"File {json_file} not found")

        if json_file.endswith(".gz"):
            with gzip.open(json_file, "rt") as f:
                json_data = json.load(f)
        else:
            with open(json_file, "r") as f:
                json_data = json.load(f)

        return cls(
            cnv=CNVRegionAnnotation(**json_data["cnv"]),
            _benign_cnv=json_data["_benign_cnv"],
            _benign_cnv_gs_inner=json_data["_benign_cnv_gs_inner"],
            _benign_cnv_gs_outer=json_data["_benign_cnv_gs_outer"],
            _regulatory=json_data["_regulatory"],
            _gnomad=json_data["_gnomad"],
            _hi_gene=json_data["_hi_gene"],
            _hi_region=json_data["_hi_region"],
            _genes=json_data["_genes"],
        )

    def get_gene_by_name(self, gene_name: str) -> dict[str, Any] | None:
        """Get 'Genes' entry with gene_name"""
        for gene in self._genes:
            if gene.get("gene_name", "") == gene_name:
                return gene
        return None

    # TODO gene types should be some enum?
    def get_genes(
        self, gene_type: str | None = None, overlap: enums.OverlapType = enums.OverlapType.all
    ) -> list[dict[str, Any]]:
        """Get 'Genes' entries whose gene_type is 'gene_type' and overlaps the CNV region in that manner"""
        if gene_type:
            genes = [gene for gene in self._genes if gene["gene_type"] == gene_type]
        else:
            genes = self._genes
        return [gene for gene in genes if self.cnv.is_overlapping(gene["start"], gene["end"], overlap)]

    def get_enhancers(self) -> list[dict[str, Any]]:
        """Get 'Regulatory' entries whose type is 'enhancer'"""
        return [region for region in self._regulatory if region["type"] == "enhancer"]

    def count_gene_types(self) -> GenesDBGeneTypesCounter:
        counter: GenesDBGeneTypesCounter = {
            "protein_coding": 0,
            "pseudogenes": 0,
            "lncrna": 0,
            "rrna": 0,
            "snrna": 0,
            "mirna": 0,
            "other": 0,
        }

        for gene in self.get_genes():
            gene_type = gene.get("gene_type", "")
            if GenesDBGeneTypes.PROTEIN_CODING in gene_type:
                counter["protein_coding"] += 1
            elif GenesDBGeneTypes.PSEUDOGENE in gene_type:
                counter["pseudogenes"] += 1
            elif GenesDBGeneTypes.LINC_RNA in gene_type:
                counter["lncrna"] += 1
            elif GenesDBGeneTypes.R_RNA in gene_type:
                counter["rrna"] += 1
            elif GenesDBGeneTypes.S_NRNA in gene_type:
                counter["snrna"] += 1
            elif GenesDBGeneTypes.MIRNA in gene_type:
                counter["mirna"] += 1
            else:
                counter["other"] += 1
        return counter

    def count_regulatory_types(self) -> RegulatoryTypesCounter:
        counter: RegulatoryTypesCounter = {
            "enhancer": 0,
            "promoter": 0,
            "open_chromatin_region": 0,
            "CTCF_binding_site": 0,
            "TF_binding_site": 0,
            "curated": 0,
            "flanking_region": 0,
            "silencer": 0,
            "transcriptional_cis_regulatory_region": 0,
            "DNase_I_hypersensitive_site": 0,
            "enhancer_blocking_element": 0,
            "TATA_box": 0,
            "other": 0,
        }
        regulatory_types = [region["type"] for region in self._regulatory]
        for reg_type in regulatory_types:
            if reg_type == RegulatoryTypes.ENHANCER:
                counter["enhancer"] += 1
            elif reg_type == RegulatoryTypes.PROMOTER:
                counter["promoter"] += 1
            elif reg_type == RegulatoryTypes.OPEN_CHROMATIN_REGION:
                counter["open_chromatin_region"] += 1
            elif reg_type == RegulatoryTypes.CTCF_BINDING_SITE:
                counter["CTCF_binding_site"] += 1
            elif reg_type == RegulatoryTypes.TF_BINDING_SITE:
                counter["TF_binding_site"] += 1
            elif reg_type == RegulatoryTypes.CURATED:
                counter["curated"] += 1
            elif reg_type == RegulatoryTypes.FLANKING_REGION:
                counter["flanking_region"] += 1
            elif reg_type == RegulatoryTypes.SILENCER:
                counter["silencer"] += 1
            elif reg_type == RegulatoryTypes.TRANSCRIPTIONAL_CIS_REGULATORY_REGION:
                counter["transcriptional_cis_regulatory_region"] += 1
            elif reg_type == RegulatoryTypes.DNASE_I_HYPERSENSITIVE_SITE:
                counter["DNase_I_hypersensitive_site"] += 1
            elif reg_type == RegulatoryTypes.ENHANCER_BLOCKING_ELEMENT:
                counter["enhancer_blocking_element"] += 1
            elif reg_type == RegulatoryTypes.TATA_BOX:
                counter["TATA_box"] += 1
            else:
                counter["other"] += 1
        return counter

    @cached_property
    def _ts_regions(self) -> list[dict[str, Any]]:
        return [
            reg
            for reg in self._hi_region
            if str(reg.get("Triplosensitivity Score", "")) in core.sufficient_HI_TS_scores
        ]

    @cached_property
    def _hi_regions(self) -> list[dict[str, Any]]:
        return [
            reg
            for reg in self._hi_region
            if str(reg.get("Haploinsufficiency Score", "")) in core.sufficient_HI_TS_scores
        ]

    @cached_property
    def _ts_genes(self) -> list[dict[str, Any]]:
        return [
            reg for reg in self._hi_gene if str(reg.get("Triplosensitivity Score", "")) in core.sufficient_HI_TS_scores
        ]

    @cached_property
    def _hi_genes(self) -> list[dict[str, Any]]:
        return [
            reg for reg in self._hi_gene if str(reg.get("Haploinsufficiency Score", "")) in core.sufficient_HI_TS_scores
        ]

    def get_ts_regions(self, overlap_type: enums.OverlapType | None) -> list[dict[str, Any]]:
        """Get all HiRegion entries with sufficient Triplosensitivity Score"""
        if not overlap_type:
            return self._ts_regions
        else:
            return [reg for reg in self._ts_regions if self.cnv.is_overlapping(reg["start"], reg["end"], overlap_type)]

    def get_hi_regions(self, overlap_type: enums.OverlapType | None) -> list[dict[str, Any]]:
        """Get all HiRegion entries with sufficient Haploinsufficiency Score"""
        if not overlap_type:
            return self._hi_regions
        else:
            return [reg for reg in self._hi_regions if self.cnv.is_overlapping(reg["start"], reg["end"], overlap_type)]

    def get_ts_genes(self, overlap_type: enums.OverlapType | None) -> list[dict[str, Any]]:
        """Get all HiGene entries with sufficient Triplosensitivity Score"""
        if not overlap_type:
            return self._ts_genes
        else:
            return [
                gene for gene in self._ts_genes if self.cnv.is_overlapping(gene["start"], gene["end"], overlap_type)
            ]

    def get_hi_genes(self, overlap_type: enums.OverlapType | None) -> list[dict[str, Any]]:
        """Get all HiGene entries with sufficient Haploinsufficiency Score"""
        if not overlap_type:
            return self._hi_genes
        else:
            return [
                gene for gene in self._hi_genes if self.cnv.is_overlapping(gene["start"], gene["end"], overlap_type)
            ]

    def get_common_variability_regions(self) -> list[CommonVariabilityRegion]:
        """Get GnomAD regions that intersect with the CNV and are in the NFE population"""
        COMMON_VARIABILITY_FREQUENCY_THRESHOLD = 0.01

        results: list[CommonVariabilityRegion] = []
        regions = [
            variability
            for variability in self._gnomad
            if variability["start"] <= self.cnv.start
            and variability["end"] >= self.cnv.end
            and variability["population"] == "nfe"
            and _normalize_type(variability["svtype"]) == self.cnv.cnv_type
        ]

        for region in regions:
            homref = region["frequencies"]["all"]["HOMREF"]["count"]
            homalt = region["frequencies"]["all"]["HOMALT"]["count"]
            het = region["frequencies"]["all"]["HET"]["count"]
            if homalt + homref + het == 0:
                continue
            freq = (homalt + het) / float(homalt + homref + het)
            if freq > COMMON_VARIABILITY_FREQUENCY_THRESHOLD:
                results.append({"population": region["population"], "frequency": freq})
        return results

    def get_benign_cnvs_gs_outer(self) -> list[dict[str, Any]]:
        """Get Benign CNV GS outer entries with no defined frequency or frequency high enough"""
        cnvs = [
            cnv
            for cnv in self._benign_cnv_gs_outer
            if cnv.get("frequency", "") == ""
            or math.isnan(cnv["frequency"])
            or float(cnv["frequency"]) >= float(core.MIN_FREQUENCY_BENIGN)
        ]

        if self.cnv.is_duplication:
            return [cnv for cnv in cnvs if self.cnv.matches_cnv_type_or_both(cnv["cnv_type"])]
        else:
            return cnvs  # TODO does not make any sense?

    def get_gene_transcript_regions(self, gene_name: str) -> list[TranscriptRegion]:
        gene_info = self.get_gene_by_name(gene_name)
        if gene_info is None:
            return []

        results: list[TranscriptRegion] = []

        is_forward = gene_info["strand"] == "+"
        transcripts = [
            t for t in gene_info.get("transcript", []) if self.cnv.get_overlap_with_region(t["start"], t["end"]) > 0
        ]

        for transcript in transcripts:
            five_inside = False
            three_inside = False
            contained = False
            cds = [
                length
                for t in transcript.get("CDS", [])
                if (length := self.cnv.get_overlap_with_region(t["start"], t["end"]))
            ]
            five_prime_utr = [
                length
                for t in transcript.get("five_prime_UTR", [])
                if (length := self.cnv.get_overlap_with_region(t["start"], t["end"]))
            ]
            three_prime_utr = [
                length
                for t in transcript.get("three_prime_UTR", [])
                if (length := self.cnv.get_overlap_with_region(t["start"], t["end"]))
            ]
            if (is_forward and self.cnv.is_position_inside(transcript["start"])) or (
                not is_forward and self.cnv.is_position_inside(transcript["end"])
            ):
                five_inside = True
            if (is_forward and self.cnv.is_position_inside(transcript["end"])) or (
                not is_forward and self.cnv.is_position_inside(transcript["start"])
            ):
                three_inside = True

            if self.cnv.is_position_inside(transcript["start"]) and self.cnv.is_position_inside(transcript["end"]):
                contained = True
            results.append(
                {
                    "length": transcript["end"] - transcript["start"],
                    "identifier": transcript["ID"],
                    "cds_overlaps": cds,
                    "five_prime_utr_overlaps": five_prime_utr,
                    "three_prime_utr_overlaps": three_prime_utr,
                    "flag_five_inside": five_inside,
                    "flag_three_inside": three_inside,
                    "flag_contained": contained,
                }
            )
        return results

    def get_high_risk_loss_genes(self) -> list[dict[str, Any]]:
        if self.cnv.is_duplication:
            raise ValueError("Only losses are considered high risk")

        high_risk_genes: list[dict[str, Any]] = []
        minimum_predicted = 2

        hi_predictors = ["HIPred", "Huang", "GHIS", "gnomeAD", "ExAC"]
        for gene in self.get_genes():
            risk_predictors: list[str] = []
            for pred in hi_predictors:
                try:
                    if gene["conservancy"][pred]["risk"]["loss"] == "high":
                        risk_predictors.append(pred)
                except KeyError:
                    pass
            if len(risk_predictors) >= minimum_predicted:
                high_risk_genes.append({"gene_name": gene["gene_name"], "risk_predictors": risk_predictors})

        return high_risk_genes
