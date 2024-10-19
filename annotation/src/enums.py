import enum


class CNVType(enum.StrEnum):
    LOSS = "loss"
    GAIN = "gain"


class Overlap(enum.Enum):
    """Possible types of overlap between two regions."""

    ANY = 0  # any overlap between regions
    SPAN_ENTIRE = 1  # the second region overlaps the entire first region
    START_ONLY = 2  # the second region overlap only the start of the second region (and not end)
    END_ONLY = 3  # the second region overlap only the end of the second region (and not start)
    START_OR_END = 4  # the second region overlap either only the end or the start of the second region
    CONTAINED_INSIDE = 5  # the second region is inside the entire first region


class RegulatoryType(enum.StrEnum):
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


class GenesDBGeneType(enum.StrEnum):
    """Types of genes as stored in the GenesDB database."""

    PROTEIN_CODING = "protein_coding"
    PSEUDOGENE = "pseudogene"  # NOTE there are any subtype of pseudogenes
    LINC_RNA = "lncRNA"
    R_RNA = "rRNA"
    S_NRNA = "snRNA"
    MIRNA = "miRNA"
