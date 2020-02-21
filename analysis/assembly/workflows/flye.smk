from enum import Enum
from pathlib import Path


class FlyeInputType(Enum):
    PACBIO = "--pacbio-corr"
    NANOPORE = "--nano-raw"


def infer_flye_input_type(technology: str) -> str:
    return FlyeInputType[technology.upper()].value


rule flye:
    input:
        reads = mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
        assembly       = outdir / "{sample}" / "flye" / "assembly.flye.{technology}.fasta",
        assembly_info  = outdir / "{sample}" / "flye" / "assembly_info.flye.{technology}.txt",
        assembly_graph = outdir / "{sample}" / "flye" / "assembly_graph.flye.{technology}.gfa",
    threads: 16
    shadow: "shallow"
    resources:
        mem_mb = lambda wildcards, attempt: 16000 * attempt,
    singularity: containers["flye"]
    params:
        genome_size = config["genome_size"],
        outdir = lambda wildcards, output: Path(output.assembly).parent,
        input_type = lambda wildcards: infer_flye_input_type(wildcards.technology),
        polishing_iterations = 1,
    shell:
        """
        flye {params.input_type} {input.reads} \
            --genome-size {params.genome_size} \
            --out-dir {params.outdir} \
            --threads {threads} \
            --iterations {params.polishing_iterations} 
        cat {params.outdir}/flye.log
        mv {params.outdir}/assembly.fasta {output.assembly}
        mv {params.outdir}/assembly_info.txt {output.assembly_info}
        mv {params.outdir}/assembly_graph.gfa {output.assembly_graph}
        """
