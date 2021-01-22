import sys

sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path

samplesheet = snakemake.params.samplesheet
qc_dir = Path(snakemake.params.qc_dir)


def generate_path(site: str, sample: str) -> Path:
    return qc_dir / f"subsampled/{site}/nanopore/{sample}/{sample}.subsampled.fastq.gz"


with open(snakemake.output.sample_idx, "w") as ostream:
    for idx, row in samplesheet.iterrows():
        site = row["site"]
        sample = row["sample"]
        path = generate_path(site, sample).resolve(strict=True)
        print(f"{sample}\t{path}", file=ostream)
