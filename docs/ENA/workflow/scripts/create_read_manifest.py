"""https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html#manifest-file"""
import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
import pandas as pd

samplesheet = pd.read_csv(snakemake.input.info, index_col="sample")
sample_acc_df = pd.read_csv(snakemake.input.samples, index_col="alias")
data_dir = Path(snakemake.params.data_dir)
reads_dir = Path(snakemake.params.reads_dir)
project = snakemake.params.project
name = snakemake.wildcards.sample
sample = sample_acc_df.at[name, "biosample"]
site = samplesheet.at[name, "site"]
tech = snakemake.wildcards.tech

if tech == "nanopore":
    platform = "OXFORD_NANOPORE"
    gridion = set(snakemake.params.gridio)
    run = samplesheet.at[name, "run"]
    instrument = "GridION" if run in gridion else "MinION"
    fastq = [reads_dir / f"{site}/{tech}/{name}/{name}.subsampled.fastq.gz"]
elif tech == "pacbio":
    platform = "PACBIO_SMRT"
    instrument = "Sequel"
    fastq = [data_dir / f"{site}/{tech}/{name}/{name}.pacbio.fastq.gz"]
elif tech == "illumina":
    platform = "ILLUMINA"
    if site == "madagascar":
        instrument = "Illumina HiSeq 2500"
        insert_size = "250"
    elif site == "south_africa":
        nextseq = set(snakemake.params.nextseq)
        instrument = "NextSeq 500" if name in nextseq else "Illumina HiSeq 2500"
        insert_size = "150" if name in nextseq else "250"
    elif site == "birmingham":
        instrument = "Illumina MiSeq"
        insert_size = "300"
    else:
        raise KeyError(f"Unrecognised site {site}")
    fastq = [
        reads_dir / f"{site}/{tech}/{name}/{name}.subsampled.R1.fastq.gz",
        reads_dir / f"{site}/{tech}/{name}/{name}.subsampled.R2.fastq.gz",
    ]
else:
    raise ValueError(f"Unknown technology {tech}")

library_selection = "unspecified"
library_source = "GENOMIC"
library_strategy = "WGS"

with open(snakemake.output.manifest, "w") as fp:
    print(f"STUDY\t{project}", file=fp)
    print(f"SAMPLE\t{sample}", file=fp)
    print(f"NAME\t{name}", file=fp)
    print(f"INSTRUMENT\t{instrument}", file=fp)
    if tech == "illumina":
        print(f"INSERT_SIZE\t{insert_size}", file=fp)

    print(f"LIBRARY_SOURCE\t{library_source}", file=fp)
    print(f"LIBRARY_SELECTION\t{library_selection}", file=fp)
    print(f"LIBRARY_STRATEGY\t{library_strategy}", file=fp)
    for p in fastq:
        print(f"FASTQ\t{p}", file=fp)
