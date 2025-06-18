import gzip
import json
import os
from dataclasses import dataclass
from typing import Any, TypedDict

from annotation.src import core, entry_schemas, enums, exceptions


class GenesDBGeneTypesCounter(TypedDict):
    protein_coding: int
    pseudogenes: int
    lncrna: int
    rrna: int
    snrna: int
    mirna: int
    other: int


class AnnotatedGenesList(TypedDict):
    morbid_genes: list[str]
    morbid_genes_urls: list[str]
    associated_with_disease: list[str]
    associated_with_disease_urls: list[str]


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
    """Representation of a transcript region"""

    length: int
    """ Length of the transcript """
    identifier: str
    """ Identifier of the transcript """
    cds_overlaps: list[int]
    """ List of lengths of CDS overlaps, in the number of bases """
    five_prime_utr_overlaps: list[int]
    """ List of lengths of prime 5' UTR overlaps, in the number of bases """
    three_prime_utr_overlaps: list[int]
    """ List of lengths of prime 3' UTR overlaps, in the number of bases """

    flag_five_inside: bool
    """ Flag indicating if the transcript start is inside the CNV region. Properly handling +/- strands"""
    flag_three_inside: bool
    """ Flag indicating if the transcript end is inside the CNV region. Properly handling +/- strands"""
    flag_contained: bool
    """ Flag indicating if the transcript is completely contained inside the CNV region """


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
    """ Chromosome of the CNV region in the format like `chr1`"""
    start: int
    """ Start position of the CNV region """
    end: int
    """ End position of the CNV region """
    cnv_type: enums.CNVType
    """ Type of the CNV region """
    length: int
    """ Length of the CNV region """
    name: str
    """ Name of the CNV region given as `{chr}_{start}_{end}_{cnv_type}`"""
    cytogenetic_position: str
    """ Cytogenetic position of the CNV region, e.g.: `q15.1`"""

    @property
    def is_duplication(self) -> bool:
        """Check if the CNV region is a duplication"""
        return self.cnv_type == enums.CNVType.GAIN

    @property
    def genomic_coord(self) -> str:
        """Genomic coordinates in the format `chr:start-end`"""
        return f"{self.chr}:{self.start}-{self.end}"

    def matches_cnv_type_or_both(self, benign_cnv_type: str) -> bool:
        """
        Check whether benign_cnv matches the CNV type or is a combination of both

        Parameters
        ----------
        benign_cnv_type : str
            Type of CNV as given by entry in benign CNV DB

        Returns
        -------
        bool
            True if benign_cnv_type matches the CNV type or is a combination of both, False otherwise
        """
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

    def is_overlapping(self, target_start: int, target_end: int, overlap: enums.Overlap) -> bool:
        """Check if the CNV region overlaps with another region in the specified manner"""
        if overlap == enums.Overlap.ANY:
            return target_start < self.end and target_end > self.start
        elif overlap == enums.Overlap.SPAN_ENTIRE:
            return target_start <= self.start and target_end >= self.end
        elif overlap == enums.Overlap.START_ONLY:
            return target_start <= self.start < target_end <= self.end
        elif overlap == enums.Overlap.END_ONLY:
            return self.start < target_start <= self.end <= target_end
        elif overlap == enums.Overlap.START_OR_END:
            return (target_start <= self.start < target_end < self.end) or (
                self.start < target_start < self.end <= target_end
            )
        elif overlap == enums.Overlap.CONTAINED_INSIDE:
            return self.start <= target_start and target_end <= self.end
        else:
            raise exceptions.UnknownOverlapTypeError(overlap, enums.Overlap.values())

    def get_overlap_with_region(self, target_start: int, target_end: int) -> int:
        """Get the number of overlapping bases with another interval."""
        return max(0, min(target_end, self.end) - max(target_start, self.start))

    def is_position_inside(self, position: int) -> bool:
        """Check if the position is inside the CNV region"""
        return self.start <= position <= self.end


@dataclass(frozen=True)
class Annotation:
    cnv: CNVRegionAnnotation
    """ See [`CNVRegionAnnotation`][src.annotation.CNVRegionAnnotation] """
    _benign_cnv: list[entry_schemas.BenignCNVDBEntry]
    _benign_cnv_gs_inner: list[entry_schemas.BenignCNVGSDBEntry]
    _benign_cnv_gs_outer: list[entry_schemas.BenignCNVGSDBEntry]
    _regulatory: list[entry_schemas.RegulatoryEntry]
    _gnomad: list[entry_schemas.GnomADEntry]
    _hi_gene: list[entry_schemas.HiGeneEntry]
    _hi_region: list[entry_schemas.HiRegionEntry]
    _genes: list[entry_schemas.GenesDBEntry]

    @classmethod
    def load_from_json(cls, json_file: str) -> "Annotation":
        """
        Load annotation from a JSON file

        Parameters
        ----------
        json_file : str
            Path to the JSON file. Can be gzipped.

        Returns
        -------
        Annotation
            Annotation object loaded from the JSON file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        """
        if not os.path.exists(json_file):
            raise FileNotFoundError(f"File {json_file} not found")

        if json_file.endswith(".gz"):
            with gzip.open(json_file, "rt") as f:
                json_data = json.load(f)
        else:
            with open(json_file, "r") as f:
                json_data = json.load(f)

        try:
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
        except KeyError as e:
            raise exceptions.InvalidJSONAnnotationError(json_file, str(e))

    def get_gene_by_name(self, gene_name: str) -> entry_schemas.GenesDBEntry | None:
        """Get 'Genes' entry with gene_name"""
        for gene in self._genes:
            if gene.get("gene_name", "") == gene_name:
                return gene
        return None

    # TODO gene types should be some enum?
    def get_genes(
        self, gene_type: str | None = None, overlap: enums.Overlap = enums.Overlap.ANY
    ) -> list[entry_schemas.GenesDBEntry]:
        """Get 'Genes' entries whose gene_type is 'gene_type' and overlaps the CNV region in that manner"""
        if gene_type:
            genes = [gene for gene in self._genes if gene["gene_type"] == gene_type]
        else:
            genes = self._genes
        return [gene for gene in genes if self.cnv.is_overlapping(gene["start"], gene["end"], overlap)]

    def get_annotated_genes(self) -> AnnotatedGenesList:
        """
        Get list of names for genes that are morbid or associated with a disease

        Returns
        -------
        AnnotatedGenesList
            Dictionary with two lists of gene names: 'morbid_genes' and 'associated_with_disease'
            and two list s of URLs: 'morbid_genes_url' and 'associated_with_disease_url'.
        """
        genes = self.get_genes()
        annot_genes: AnnotatedGenesList = {
            "morbid_genes": [],
            "morbid_genes_urls": [],
            "associated_with_disease": [],
            "associated_with_disease_urls": [],
        }

        for gene in genes:
            if "AnnotSV" in gene:
                if gene["AnnotSV"].get("omim_morbid_gene", "") == "yes":
                    annot_genes["morbid_genes"].append(gene["gene_name"])
                    annot_genes["morbid_genes_urls"].append(gene["external"]["OMIM"]["url"])
                if "omim_phenotype" in gene["AnnotSV"]:
                    annot_genes["associated_with_disease"].append(gene["gene_name"])
                    annot_genes["associated_with_disease_urls"].append(gene["external"]["OMIM"]["url"])
        return annot_genes

    def get_triplosensitivity_regions(
        self, overlap_type: enums.Overlap, valid_scores: list[int]
    ) -> list[entry_schemas.HiRegionEntry]:
        """
        Get all HiRegion entries with sufficient Triplosensitivity Score

        Parameters
        ----------
        overlap_type : enums.Overlap
            Type of overlap to consider.
        valid_scores : list[int]
            List of valid scores to consider. Entry scores are coerced to integers.

        Returns
        -------
        list[entry_schemas.HiRegionEntry]
            List of HiRegion entries with sufficient Triplosensitivity Score
        """
        return [
            r
            for r in self._hi_region
            if "Triplosensitivity Score" in r
            and r["Triplosensitivity Score"] not in core.INVALID_SCORES
            and int(float(r["Triplosensitivity Score"])) in valid_scores
            and self.cnv.is_overlapping(r["start"], r["end"], overlap_type)
        ]

    def get_haploinsufficient_regions(
        self, overlap_type: enums.Overlap, valid_scores: list[int]
    ) -> list[entry_schemas.HiRegionEntry]:
        """
        Get all HiRegion entries with sufficient Haploinsufficiency Score

        Parameters
        ----------
        overlap_type : enums.Overlap
            Type of overlap to consider.
        valid_scores : list[int]
            List of valid scores to consider. Entry scores are coerced to integers.

        Returns
        -------
        list[entry_schemas.HiRegionEntry]
            List of HiRegion entries with sufficient Haploinsufficiency Score
        """
        return [
            r
            for r in self._hi_region
            if "Haploinsufficiency Score" in r
            and r["Haploinsufficiency Score"] not in core.INVALID_SCORES
            and int(float(r["Haploinsufficiency Score"])) in valid_scores
            and self.cnv.is_overlapping(r["start"], r["end"], overlap_type)
        ]

    def get_triplosensitivity_genes(
        self, overlap_type: enums.Overlap, valid_scores: list[int]
    ) -> list[entry_schemas.HiGeneEntry]:
        """
        Get all HiGene entries with sufficient Triplosensitivity Score

        Parameters
        ----------
        overlap_type : enums.Overlap
            Type of overlap to consider.
        valid_scores : list[int]
            List of valid scores to consider. Entry scores are coerced to integers.

        Returns
        -------
        list[entry_schemas.HiGeneEntry]
            List of HiGene entries with sufficient Triplosensitivity Score
        """
        return [
            g
            for g in self._hi_gene
            if "Triplosensitivity Score" in g
            and g["Triplosensitivity Score"] not in core.INVALID_SCORES
            and int(float(g["Triplosensitivity Score"])) in valid_scores
            and self.cnv.is_overlapping(g["start"], g["end"], overlap_type)
        ]

    def get_triplosensitivity_gene_names(self, overlap_type: enums.Overlap, valid_scores: list[int]) -> list[str]:
        """Get names of genes with sufficient Triplosensitivity Score"""
        return [gene["Gene Symbol"] for gene in self.get_triplosensitivity_genes(overlap_type, valid_scores)]

    def get_haploinsufficient_genes(
        self, overlap_type: enums.Overlap, valid_scores: list[int]
    ) -> list[entry_schemas.HiGeneEntry]:
        """
        Get all HiGene entries with sufficient Haploinsufficiency Score

        Parameters
        ----------
        overlap_type : enums.Overlap
            Type of overlap to consider.
        valid_scores : list[int]
            List of valid scores to consider. Entry scores are coerced to integers.

        Returns
        -------
        list[entry_schemas.HiGeneEntry]
            List of HiGene entries with sufficient Haploinsufficiency Score
        """
        return [
            g
            for g in self._hi_gene
            if "Haploinsufficiency Score" in g
            and g["Haploinsufficiency Score"] not in core.INVALID_SCORES
            and int(float(g["Haploinsufficiency Score"])) in valid_scores
            and self.cnv.is_overlapping(g["start"], g["end"], overlap_type)
        ]

    def get_haploinsufficient_gene_names(self, overlap_type: enums.Overlap, valid_scores: list[int]) -> list[str]:
        """Get names of genes with sufficient Haploinsufficiency Score"""
        return [gene["Gene Symbol"] for gene in self.get_haploinsufficient_genes(overlap_type, valid_scores)]

    def get_hi_or_ts_genes_url(self, hi_or_genes_list: list[str]) -> list[str]:

        "Get URLs for genes that are in list of Haploinsufficient or Triplosensitive genes"
        data = self._genes

        hi_or_ts_gene_mane_url = []

        for hi_or_ts_gene_mane in hi_or_genes_list:

            for i in range(len(data)):

                if data[i]['gene_name'] == hi_or_ts_gene_mane:
                    if 'external' in data[i] and 'OMIM' in data[i]['external']:
                        omim_url = data[i]['external']['OMIM']['url']
                        hi_or_ts_gene_mane_url.append(omim_url)
                    else:
                        print("No OMIM information available for this gene.")
                    break

        return hi_or_ts_gene_mane_url

    def get_common_variability_regions(self) -> list[CommonVariabilityRegion]:
        """
        Get GnomAD regions that intersect with the CNV and are in the specified populations

        Returns
        -------
        list[CommonVariabilityRegion]
            List of common variability regions that intersect with the CNV

        See Also
        --------
        core.POPULATION_FOR_COMMON_VARIABILITY
            List of populations for common variability
        core.COMMON_VARIABILITY_FREQUENCY_THRESHOLD
            Minimum frequency for a common variability region to be considered
        """

        results: list[CommonVariabilityRegion] = []
        regions = [
            variability
            for variability in self._gnomad
            if variability["start"] <= self.cnv.start
            and variability["end"] >= self.cnv.end
            and variability["population"] in core.POPULATION_FOR_COMMON_VARIABILITY
            and _normalize_type(variability["svtype"]) == self.cnv.cnv_type
        ]

        for region in regions:
            homref = region["frequencies"]["all"]["HOMREF"]["count"]
            homalt = region["frequencies"]["all"]["HOMALT"]["count"]
            het = region["frequencies"]["all"]["HET"]["count"]
            if homalt + homref + het == 0:
                continue
            freq = (homalt + het) / float(homalt + homref + het)
            if freq > core.COMMON_VARIABILITY_FREQUENCY_THRESHOLD:
                results.append({"population": region["population"], "frequency": freq})
        return results

    def get_benign_cnvs_gs_outer(self, frequency_threshold: float) -> list[entry_schemas.BenignCNVGSDBEntry]:
        """
        Get Benign CNV GS outer entries with no defined frequency or higher than the threshold

        Parameters
        ----------
        frequency_threshold : float
            Minimum frequency for a benign CNV GS outer entry to be considered

        Returns
        -------
        list[dict[str, Any]]
            List of benign CNV GS outer entries

        See Also
        --------
        core.MIN_FREQUENCY_BENIGN
            Minimum frequency for a benign CNV GS outer entry to be considered
        """
        cnvs = [
            cnv
            for cnv in self._benign_cnv_gs_outer
            if "frequency" not in cnv
            or (isinstance(cnv["frequency"], str) and cnv["frequency"].lower() == "nan")
            or float(cnv["frequency"]) >= frequency_threshold
        ]  # TODO check possible types for frequency in schema

        if self.cnv.is_duplication:
            return [cnv for cnv in cnvs if self.cnv.matches_cnv_type_or_both(cnv["cnv_type"])]
        else:
            return cnvs  # TODO does not make any sense?

    def get_gene_transcript_regions(self, gene_name: str) -> list[TranscriptRegion]:
        """
        Get transcript regions for a gene

        Parameters
        ----------
        gene_name : str
            Name of the gene

        Returns
        -------
        list[TranscriptRegion]
            List of transcript regions for the gene
        """
        gene_info = self.get_gene_by_name(gene_name)
        if gene_info is None:
            return []

        results: list[TranscriptRegion] = []

        is_forward = gene_info["strand"] == "+"
        transcripts: list[entry_schemas.TranscriptEntry] = [
            t for t in gene_info.get("transcript", {}) if self.cnv.get_overlap_with_region(t["start"], t["end"]) > 0
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
        """
        Get genes with high risk loss predictors

        Raises
        ------
        HighRiskForDuplicationError
            If the CNV region is a duplication, not a loss

        Returns
        -------
        list[dict[str, Any]]
            List of genes with high risk loss predictors

        See Also
        --------
        core.HIGH_RISK_PREDICTORS
            List of high risk predictors
        core.HIGH_RISK_PREDICTORS_COUNT_THRESHOLD
            Minimum number of high risk predictors for a gene to be considered high risk
        """
        if self.cnv.is_duplication:
            raise exceptions.HighRiskForDuplicationError()

        high_risk_genes: list[dict[str, Any]] = []

        for gene in self.get_genes():
            risk_predictors: list[str] = []
            for pred in core.HIGH_RISK_PREDICTORS:
                try:
                    if gene["conservancy"][pred]["risk"]["loss"] == "high":
                        risk_predictors.append(pred)
                except KeyError:
                    pass
            if len(risk_predictors) >= core.HIGH_RISK_PREDICTORS_COUNT_THRESHOLD:
                high_risk_genes.append({"gene_name": gene["gene_name"], "risk_predictors": risk_predictors})

        return high_risk_genes

    def count_gene_types(self) -> GenesDBGeneTypesCounter:
        """
        Count the number of genes of each type in GenesDB

        Returns
        -------
        GenesDBGeneTypesCounter
            A dictionary with the count of each gene type.
        """
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
            if enums.GenesDBGeneType.PROTEIN_CODING in gene_type:
                counter["protein_coding"] += 1
            # NOTE there are any subtype of pseudogenes like 'unprocessed_pseudogene'  or 'IG_V_pseudogene'
            elif enums.GenesDBGeneType.PSEUDOGENE in gene_type:
                counter["pseudogenes"] += 1
            elif enums.GenesDBGeneType.LINC_RNA in gene_type:
                counter["lncrna"] += 1
            elif enums.GenesDBGeneType.R_RNA in gene_type:
                counter["rrna"] += 1
            elif enums.GenesDBGeneType.S_NRNA in gene_type:
                counter["snrna"] += 1
            elif enums.GenesDBGeneType.MIRNA in gene_type:
                counter["mirna"] += 1
            else:
                counter["other"] += 1
        return counter

    def count_regulatory_types(self) -> RegulatoryTypesCounter:
        """
        Count the number of regulatory element types in the database

        Returns
        -------
        RegulatoryTypesCounter
            A dictionary with the count of each regulatory type.

        Notes
        -----
            The regulatory database entries are looped over and their "type" element is checked for presence in any of existing regulatory types.
            Unknown regulatory types are counted as 'other'. To get a count of all regulatory types, sum the values.
        """

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
            if reg_type == enums.RegulatoryType.ENHANCER:
                counter["enhancer"] += 1
            elif reg_type == enums.RegulatoryType.PROMOTER:
                counter["promoter"] += 1
            elif reg_type == enums.RegulatoryType.OPEN_CHROMATIN_REGION:
                counter["open_chromatin_region"] += 1
            elif reg_type == enums.RegulatoryType.CTCF_BINDING_SITE:
                counter["CTCF_binding_site"] += 1
            elif reg_type == enums.RegulatoryType.TF_BINDING_SITE:
                counter["TF_binding_site"] += 1
            elif reg_type == enums.RegulatoryType.CURATED:
                counter["curated"] += 1
            elif reg_type == enums.RegulatoryType.FLANKING_REGION:
                counter["flanking_region"] += 1
            elif reg_type == enums.RegulatoryType.SILENCER:
                counter["silencer"] += 1
            elif reg_type == enums.RegulatoryType.TRANSCRIPTIONAL_CIS_REGULATORY_REGION:
                counter["transcriptional_cis_regulatory_region"] += 1
            elif reg_type == enums.RegulatoryType.DNASE_I_HYPERSENSITIVE_SITE:
                counter["DNase_I_hypersensitive_site"] += 1
            elif reg_type == enums.RegulatoryType.ENHANCER_BLOCKING_ELEMENT:
                counter["enhancer_blocking_element"] += 1
            elif reg_type == enums.RegulatoryType.TATA_BOX:
                counter["TATA_box"] += 1
            else:
                counter["other"] += 1
        return counter
