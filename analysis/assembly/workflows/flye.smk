from enum import Enum
from pathlib import Path


class FlyeInputType(Enum):
    PACBIO = "--pacbio-hifi"
    NANOPORE = "--nano-raw"


def infer_flye_input_type(technology: str) -> str:
    return FlyeInputType[technology.upper()].value


rule flye:
    input:
        reads=mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
        assembly=outdir / "{sample}" / "flye" / "{technology}" / "assembly.flye.{technology}.fasta",
        assembly_info=outdir / "{sample}" / "flye" / "{technology}" / "assembly_info.flye.{technology}.txt",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
    singularity: containers["flye"]
    params:
        genome_size=config["genome_size"],
        outdir=lambda wildcards, output: Path(output.assembly).parent,
        input_type=lambda wildcards: infer_flye_input_type(wildcards.technology),
        polishing_iterations=1,
        original_assembly_info_name="assembly_info.txt",
        original_assembly_name="assembly.fasta",
    shell:
        """
        flye {params.input_type} {input.reads} \
            --genome-size {params.genome_size} \
            --out-dir {params.outdir} \
            --threads {threads} \
            --iterations {params.polishing_iterations} 
        sleep 5
        cp {params.outdir}/{params.original_assembly_name} {output.assembly}
        cp {params.outdir}/{params.original_assembly_info_name} {output.assembly_info}
        """
