from enum import Enum

class FlyeInputType(Enum):
    PACBIO = "--pacbio-corr"
    NANOPORE = "--nano-raw"


def infer_flye_input_type(technology: str) -> str:
    return FlyeInputType[technology.upper()].value


rule flye:
    input:
        reads = mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
        assembly      = outdir / "{sample}" / "flye" / "{technology}" / "assembly.fasta",
        assembly_info = outdir / "{sample}" / "flye" / "{technology}" / "assembly_info.txt",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 16000 * attempt,
    singularity: "docker://quay.io/biocontainers/flye:2.6--py37he513fc3_0"
    params:
        genome_size = config["genome_size"],
        input_type = lambda wildcards: infer_flye_input_type(wildcards.technology),
        polishing_iterations = 1,
    shell:
        """
        outdir=$(dirname {output.assembly})
        flye {params.input_type} {input.reads} \
            --genome-size {params.genome_size} \
            --out-dir $outdir \
            --threads {threads} \
            --iterations {params.polishing_iterations} 
        """
