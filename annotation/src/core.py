import os

SRC_DIR = os.path.abspath(os.path.dirname(__file__))

CYTOBANDS_FILEPATH = os.path.join(SRC_DIR, "cytobands", "cytoBand_hg38.txt")


REGION_REGEX = r"(chr[\dXY]+):(\d+)-(\d+)/(\w+)"
JSON_INDENT_LEVEL = 2
ENCODING = "utf-8"

INVALID_SCORES = ("", "Not yet evaluated", "nan")

# TODO rework to use as params
HIGH_RISK_PREDICTORS = ["HIPred", "Huang", "GHIS", "gnomeAD", "ExAC"]
HIGH_RISK_PREDICTORS_COUNT_THRESHOLD = 2
COMMON_VARIABILITY_FREQUENCY_THRESHOLD = 0.01
POPULATION_FOR_COMMON_VARIABILITY = ["nfe"]
