import sys
from typing import Tuple

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import pysam
from pathlib import Path


def asm_stats(file: Path) -> Tuple[int, int]:
    n_seqs = 0
    bases = 0
    with pysam.FastxFile(file) as fp:
        for entry in fp:
            n_seqs += 1
            bases += len(entry.sequence)
    return n_seqs, bases


def coverage(fastq: Path, asm_length: int) -> int:
    bases_in_reads = 0
    with pysam.FastxFile(fastq) as fq:
        for entry in fq:
            bases_in_reads += len(entry.sequence)

    return int(bases_in_reads / asm_length)


def create_chromosome_list(asm_file: Path, outpath: Path):
    with pysam.FastxFile(asm_file) as fa:
        entry = next(fa)
        obj_name = entry.name
        chrom_name = "1"
        chrom_type = "Circular-Chromosome"

    with open(outpath, "w") as ofp:
        print(f"{obj_name}\t{chrom_name}\t{chrom_type}", file=ofp)


asm_path = Path(snakemake.input.asm)
data_dir = Path(snakemake.params.data_dir)
alias = snakemake.wildcards.sample
study = snakemake.params["project"]
sample_acc_df = pd.read_csv(snakemake.input.samples, index_col="alias")
run_acc_df = pd.read_csv(snakemake.input.runs, index_col=["alias", "tech"])
sample = sample_acc_df.at[alias, "biosample"]
run = run_acc_df.at[(alias, "pacbio"), "run"]

n_contigs, n_bases = asm_stats(asm_path)

if n_contigs == 1:
    chrom_list_file = asm_path.parent / f"{alias}.chrom_list.tsv"
    create_chromosome_list(asm_path, chrom_list_file)

fastq = data_dir / f"madagascar/pacbio/{alias}/{alias}.pacbio.fastq.gz"
covg = coverage(fastq, n_bases)

with open(snakemake.output.manifest, "w") as fp:
    # STUDY: Study accession - mandatory
    print(f"STUDY\t{study}", file=fp)
    # SAMPLE: Sample accession - mandatory
    print(f"SAMPLE\t{sample}", file=fp)
    # ASSEMBLYNAME: Unique assembly name, user-provided - mandatory
    print(f"ASSEMBLYNAME\t{alias}", file=fp)
    # ASSEMBLY_TYPE: ‘clone or isolate’ - mandatory
    print("ASSEMBLY_TYPE\tisolate", file=fp)
    # COVERAGE: The estimated depth of sequencing coverage - mandatory
    print(f"COVERAGE\t{covg}", file=fp)
    # PROGRAM: The assembly program - mandatory
    print("PROGRAM\tflye", file=fp)
    # PLATFORM: The sequencing platform, or comma-separated list of platforms - mandatory
    print("PLATFORM\tPACBIO_SMRT", file=fp)
    # MINGAPLENGTH: Minimum length of consecutive Ns to be considered a gap - optional
    # MOLECULETYPE: ‘genomic DNA’, ‘genomic RNA’ or ‘viral cRNA’ - optional
    print("MOLECULETYPE\tgenomic DNA", file=fp)
    # DESCRIPTION: Free text description of the genome assembly - optional
    # RUN_REF: Comma separated list of run accession(s) - optional
    print(f"RUN_REF\t{run}", file=fp)
    # FASTA: sequences in fasta format
    print(f"FASTA\t{asm_path.resolve()}", file=fp)
    # CHROMOSOME_LIST: list of chromosomes
    if n_contigs == 1:
        print(f"CHROMOSOME_LIST\t{chrom_list_file.resolve()}")
