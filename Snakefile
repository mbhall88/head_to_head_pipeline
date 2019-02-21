import re
from pathlib import Path
from typing import List

import pandas as pd
from snakemake.utils import min_version, validate

min_version("5.1.0")


# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"
validate(config, "analysis/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"]).set_index(
    ["region", "nanopore_run_id", "sample_id"], drop=False
)
validate(samples, "analysis/schemas/samples.schema.yaml")
# ======================================================
# Functions and Classes
# ======================================================


# ======================================================
# Global variables
# ======================================================
RULES_DIR = Path("analysis/rules")

files = []
for index, row in samples.iterrows():
    run_id = row["nanopore_run_id"]
    if pd.isnull(run_id):
        continue
    region = row["region"]
    sample_id = row["sample_id"]
    files.extend([
        "analysis/{region}/nanopore/{run}/report/{sample}.report.html".format(
            region=region, run=run_id, sample=sample_id
        )])


# ======================================================
# Rules
# ======================================================
rule all:
    input:
        files

# the snakemake files that run the different parts of the pipeline
include: str(RULES_DIR.joinpath("basecall.smk"))
include: str(RULES_DIR.joinpath("demultiplex.smk"))
include: str(RULES_DIR.joinpath("trim.smk"))
include: str(RULES_DIR.joinpath("remove_contamination.smk"))
include: str(RULES_DIR.joinpath("mykrobe.smk"))
include: str(RULES_DIR.joinpath("plotting.smk"))
include: str(RULES_DIR.joinpath("report.smk"))
