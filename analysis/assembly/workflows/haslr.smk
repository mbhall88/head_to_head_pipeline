from pathlib import Path
from enum import Enum


class HaslrInputType(Enum):
    PACBIO = "corrected"
    NANOPORE = "nanopore"


def infer_haslr_input_type(technology: str) -> str:
    return HaslrInputType[technology.upper()].value


rule haslr:
    input:
        illumina1=rules.trim_illumina.output.forward_paired,
        illumina2=rules.trim_illumina.output.reverse_paired,
        long_reads=mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
        assembly=outdir / "{sample}" / "haslr" / "{technology}" / "assembly.haslr.{technology}.fasta",
        assembly_graph=outdir / "{sample}" / "haslr" / "{technology}" / "assembly.haslr.{technology}.gfa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
    params:
        long_read_covg=0,  # use all reads
        long_read_type=lambda wildcards: infer_haslr_input_type(wildcards.technology),
        genome_size=config["genome_size"],
        outdir=lambda wildcards, output: Path(output.assembly).parent,
    singularity: containers["haslr"]
    shell:
        """
        haslr.py --out {params.outdir} \
            --genome {params.genome_size} \
            --long {input.long_reads} \
            --type {params.long_read_type} \
            --short {input.illumina1} {input.illumina2} \
            --threads {threads} \
            --cov-lr {params.long_read_covg} 
        sleep 5
        cp {params.outdir}/asm_contigs*/asm.final.fa {output.assembly}
        cp {params.outdir}/asm_contigs*/backbone.06.smallbubble.gfa {output.assembly_graph}
        """
