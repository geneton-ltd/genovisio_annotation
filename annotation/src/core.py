import os

SRC_DIR = os.path.abspath(os.path.dirname(__file__))

CYTOBANDS_FILEPATH = os.path.join(SRC_DIR, "cytobands", "cytoBand_hg38.txt")


REGION_REGEX = r"(chr[\dXY]+):(\d+)-(\d+)/(\w+)"
JSON_INDENT_LEVEL = 2
ENCODING = "utf-8"


MIN_FREQUENCY_BENIGN = 0.005

SUFFICIENT_HI_SCORES = [1, 2, 3, 30]
SUFFICIENT_TS_SCORES = ["1", "2", "3", "30"]

INVALID_HI_REGIONS_VALUES = (40, 0, "", "nan")
INVALID_TS_REGIONS_VALUES = ["40", "0", "Not yet evaluated"]

POPULATION_FOR_COMMON_VARIABILITY = ["nfe"]

HIGH_RISK_PREDICTORS = ["HIPred", "Huang", "GHIS", "gnomeAD", "ExAC"]
HIGH_RISK_PREDICTORS_COUNT_THRESHOLD = 2
COMMON_VARIABILITY_FREQUENCY_THRESHOLD = 0.01
