import re
import pandas as pd
from typing import List
from pathlib import Path
from snakemake.utils import min_version, validate


# Snakemake version when Singularity support was added
min_version("5.1.0")


# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"
validate(config, "analysis/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"]).set_index(
    ["region", "nanopore_run_id", "sample_id"], drop=False)
validate(samples, "analysis/schemas/samples.schema.yaml")


empty = ("madagscar_tb_aug_4", "madagscar_tb_aug_3", "madagascar_tb_mdr_control_7",
         "madagscar_tb_aug_2", "madagscar_tb_aug_1", "madagscar_tb_aug_6", "madagscar_tb_aug_5")
samples = samples[~samples.nanopore_run_id.isin(empty)]

# ======================================================
# Functions and Classes
# ======================================================


# ======================================================
# Global variables
# ======================================================
RULES_DIR = Path("analysis/rules")

regions = set(samples["region"])
runs = set(samples.dropna(subset=['nanopore_run_id'])["nanopore_run_id"])

basecall_fastq = set()
for index, row in samples.iterrows():
    run_id = row['nanopore_run_id']
    if pd.isnull(run_id):
        continue
    basecall_fastq.add(
        "analysis/{region}/nanopore/{run}/basecalled.fastq".format(region=row["region"], run=run_id))
# ======================================================
# Rules
# ======================================================
rule all:
    input:
        basecall_fastq

# the snakemake files that run the different parts of the pipeline
include: str(RULES_DIR.joinpath("basecall.smk"))
