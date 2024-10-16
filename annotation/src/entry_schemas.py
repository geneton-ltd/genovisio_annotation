from typing import TypedDict

HiGeneEntry = TypedDict(
    "HiGeneEntry",
    {
        "_id": str,
        "Gene ID": int,
        "cytoBand": str,
        "Genomic Location": str,
        "Haploinsufficiency Score": int,
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

GnomADFrequencyEntry = TypedDict(
    "GnomADFrequencyEntry",
    {
        "count": int,
        "frequency": float,
    },
)

GnomADGenderFrequencyEntry = TypedDict(
    "GnomADGenderFrequencyEntry",
    {
        "HET": GnomADFrequencyEntry,
        "HOMALT": GnomADFrequencyEntry,
        "HOMREF": GnomADFrequencyEntry,
    },
)

GnomADFrequenciesEntry = TypedDict(
    "GnomADFrequenciesEntry",
    {
        "male": GnomADGenderFrequencyEntry,
        "female": GnomADGenderFrequencyEntry,
        "all": GnomADGenderFrequencyEntry,
    },
)

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
