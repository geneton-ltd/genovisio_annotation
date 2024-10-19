from typing import Any, Literal, NotRequired, TypedDict

from annotation.src import enums

HiGeneEntry = TypedDict(
    "HiGeneEntry",
    {
        "_id": str,
        "Gene ID": int,
        "cytoBand": str,
        "Genomic Location": str,
        "Haploinsufficiency Score": float | str,
        "Haploinsufficiency Description": str,
        "Triplosensitivity Score": str,
        "Triplosensitivity Description": str,
        "Date Last Evaluated": str,
        "Haploinsufficiency Disease ID": str,
        "start": int,
        "end": int,
        "chromosome": str,
        "Gene Symbol": str,
    },
)


class GnomADFrequencyEntry(TypedDict):
    count: int
    frequency: float


GnomADGenderFrequencyEntry = TypedDict(
    "GnomADGenderFrequencyEntry",
    {
        "HET": GnomADFrequencyEntry,
        "HOMALT": GnomADFrequencyEntry,
        "HOMREF": GnomADFrequencyEntry,
    },
)


class GnomADFrequenciesEntry(TypedDict):
    male: GnomADGenderFrequencyEntry
    female: GnomADGenderFrequencyEntry
    all: GnomADGenderFrequencyEntry


GnomADEntry = TypedDict(
    "GnomADEntry",
    {
        "_id": str,
        "chromosome": str,
        "start": int,
        "end": int,
        "svtype": str,
        "population": str,
        "count": int,
        "frequencies": GnomADFrequenciesEntry,
    },
)


class BenignDBGeneEntry(TypedDict):
    _id: str
    chromosome: str
    start: int
    end: int
    source: str
    gene_name: str
    ID: str


class BenignCNVDBEntry(TypedDict):
    _id: str
    chromosome: str
    start: int
    end: int
    variantaccession: str
    varianttype: str
    variantsubtype: str
    reference: str
    pubmedid: int
    method: str
    supportingvariants: str
    mergedorsample: str
    samplesize: int
    observedgains: float
    observedlosses: float
    cnv_type: str
    genes: list[BenignDBGeneEntry]


class BenignCNVGSDBEntry(TypedDict):
    _id: str
    chromosome: str
    start: int
    end: int
    cnv_type: str  # TODO possible types?
    variantaccession: str
    length: int
    reference: str
    frequency: NotRequired[float | str]  # TODO can be missing or not? in Marcnv there is a check for that
    number_of_unique_samples_tested: int
    genes: list[BenignDBGeneEntry]


class TranscriptSequenceEntry(TypedDict):
    start: int
    end: int
    strand: Literal["+", "-"]
    phase: str
    ID: str
    Parent: str
    gene_id: str
    transcript_id: str
    gene_type: str
    gene_name: str
    transcript_type: str
    transcript_name: str
    exon_number: str
    exon_id: str
    level: str
    tag: str
    havana_gene: str
    havana_transcript: str


class TranscriptEntry(TypedDict):
    start: int
    end: int
    strand: Literal["+", "-"]
    phase: str
    ID: str
    Parent: str
    gene_id: str
    transcript_id: str
    gene_type: str
    gene_name: str
    transcript_type: str
    transcript_name: str
    level: str
    tag: str
    havana_gene: str
    havana_transcript: str
    exon: NotRequired[list[TranscriptSequenceEntry]]
    five_prime_UTR: NotRequired[list[TranscriptSequenceEntry]]
    CDS: NotRequired[list[TranscriptSequenceEntry]]
    three_prime_UTR: NotRequired[list[TranscriptSequenceEntry]]
    stop_codon: NotRequired[list[TranscriptSequenceEntry]]
    start_codon: NotRequired[list[TranscriptSequenceEntry]]


class GenesDBAnnotSVEntry(TypedDict):
    omim_morbid_gene: str
    omim_morbid_candidate: NotRequired[str]
    omim_phenotype: list[str]


class GenesDBEntry(TypedDict):
    _id: str
    chromosome: str
    source: str
    start: int
    end: int
    strand: Literal["+", "-"]
    phase: str
    ID: str
    gene_id: str
    gene_type: enums.GenesDBGeneType  # TODO but there are more types ...
    gene_name: str
    level: str
    havana_gene: str
    name: str
    full_name: str
    alternative_names: list[Any]  # TODO IDK if this is correct
    external: dict[str, Any]  # TODO TBS later
    pathways: list[Any]  # TODO TBS later
    function: dict[str, Any]  # TODO TBS later
    expression: dict[str, Any]  # TODO TBS later
    conservancy: dict[str, Any]  # TODO TBS later
    phenotype: dict[str, Any]  # TODO TBS later
    AnnotSV: NotRequired[GenesDBAnnotSVEntry]
    transcript: list[TranscriptEntry]


class RegulatoryEntry(TypedDict):
    _id: str
    id: str
    source: str
    type: enums.RegulatoryType
    chromosome: str
    start: int
    end: int
    ensembl_data: dict[str, Any]  # TODO TBS later
    annotations: dict[str, Any]  # TODO TBS later
    dbxrefs: list[Any]  # TODO TBS later
    features: list[dict[str, Any]]  # TODO TBS later


HiRegionEntry = TypedDict(
    "HiRegionEntry",
    {
        "_id": str,
        "ISCA Region Name": str,
        "chromosome": str,
        "start": int,
        "end": int,
        "cytoBand": str,
        "Genomic Location": str,
        "Haploinsufficiency Score": float | str,
        "Haploinsufficiency Description": str,
        "Haploinsufficiency PMID1": NotRequired[str],  # TODO what is this? there can be multiple
        "Triplosensitivity Score": str,
        "Triplosensitivity Description": str,
        "Triplosensitivity PMID1": NotRequired[str],  # TODO what is this? there can be multiple
        "Date Last Evaluated": str,
    },
)
