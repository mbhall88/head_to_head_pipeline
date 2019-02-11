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
validate(samples, "analysis/schemas/samples.schema.yaml")


empty = ("madagscar_tb_aug_4", "madagscar_tb_aug_3", "madagascar_tb_mdr_control_7", "madagscar_tb_aug_2", "madagscar_tb_aug_1", "madagscar_tb_aug_6", "madagscar_tb_aug_5")
samples = samples[~samples.nanopore_run_id.isin(empty)]
print(samples)

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

regions = set(samples["region"])
runs = set(samples.dropna(subset=['nanopore_run_id'])["nanopore_run_id"])
#======================================================
# Rules
#======================================================
rule all:
    input:
        expand("analysis/{region}/nanopore/{run}/basecalled.fastq", region=regions, run=runs)

# the snakemake files that run the different parts of the pipeline
include: str(RULES_DIR.joinpath("basecall.smk"))

