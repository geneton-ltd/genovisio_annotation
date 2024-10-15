import os

SRC_DIR = os.path.abspath(os.path.dirname(__file__))

CYTOBANDS_FILEPATH = os.path.join(SRC_DIR, "cytobands", "cytoBand_hg38.txt")


REGION_REGEX = r"(chr[\dXY]+):(\d+)-(\d+)/(\w+)"
JSON_INDENT_LEVEL = 2
ENCODING = "utf-8"

MIN_FREQUENCY_BENIGN = 0.005

sufficient_HI_TS_scores = ["3"]
