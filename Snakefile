import re
import pandas as pd
from typing import List
from pathlib import Path
from snakemake.utils import min_version, validate


# Snakemake version when Singularity support was added
min_version("5.1.0")


#======================================================
# Config files
#======================================================
configfile: "config.yaml"
validate(config, "analysis/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"]).set_index(["region", "nanopore_run_id", "sample_id"], drop=False)
samples.sort_index(inplace=True, level=2)
validate(samples, "analysis/schemas/samples.schema.yaml")

print(set(samples.xs(('madagascar', 'madagascar_tb_controls_4')).loc[ : , 'flowcell_version']))

#======================================================
# Functions and Classes
#======================================================
class InvalidBarcode(Exception):
    __module__ = Exception.__module__


def barcode_parser(barcodes_string: str) -> List[str]:
    """Parses the barcodes string and ensures they follow correct format"""
    msg = "Barcode must be of the form BC01. That is, BC followed by 2 digits."
    regex = r"\bBC\d{2}\b"
    barcodes = barcodes_string.split()
    for barcode in barcodes:
        if not (len(barcode) == 4 and re.match(regex, barcode)):
            raise InvalidBarcode(barcode + "\n" + msg)
    return barcodes


#======================================================
# Global variables
#======================================================
RULES_DIR = Path("analysis/rules")

#======================================================
# Rules
#======================================================
rule all:
    input:

# the snakemake files that run the different parts of the pipeline
include: str(RULES_DIR.joinpath("basecall.smk"))

