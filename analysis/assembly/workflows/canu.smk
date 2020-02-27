from pathlib import Path
from enum import Enum


class CanuInputType(Enum):
    PACBIO = "-pacbio-hifi"
    NANOPORE = "-nanopore-raw"


def infer_canu_input_type(technology: str) -> str:
    return CanuInputType[technology.upper()].value


rule canu:
    input:
        reads=mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
        assembly=outdir / "{sample}" / "canu" / "{technology}" / "assembly.canu.with_bubbles.{technology}.fasta",
        report=outdir / "{sample}" / "canu" / "{technology}" / "assembly.canu.with_bubbles.{technology}.report",
        assembly_graph=outdir / "{sample}" / "canu" / "{technology}" / "assembly.canu.with_bubbles.{technology}.gfa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt
    params:
        min_read_length=1000,
        min_overlap_length=500,
        genome_size=config["genome_size"],
        input_type=lambda wildcards: infer_canu_input_type(wildcards.technology),
        mem_gb=lambda wildcards, resources: int(resources.mem_mb / 1000),
        outprefix=lambda wildcards, output: Path(output.assembly).with_suffix("").name,
        outdir=lambda wildcards, output: Path(output.assembly).parent,
        extra="corMaxEvidenceErate=0.15",
        # for GC rich, repetitive genomes from docs https://canu.readthedocs.io/en/latest/faq.html#my-genome-is-at-or-gc-rich-do-i-need-to-adjust-parameters-what-about-highly-repetitive-genomes
    singularity: containers["canu"]
    shell:
        """
        canu -p {params.outprefix} \
            -d {params.outdir} \
            minReadLength={params.min_read_length} \
            minOverlapLength={params.min_overlap_length} \
            genomeSize={params.genome_size} \
            minThreads={threads} maxThreads={threads} \
            maxMemory={params.mem_gb} \
            {params.input_type} {input.reads} \
            {params.extra}
        sleep 5
        cp {params.outdir}/{params.outprefix}.contigs.fasta {output.assembly}
        cp {params.outdir}/{params.outprefix}.unitigs.gfa {output.assembly_graph}
        """

rule remove_bubbles_canu:
    input:
        assembly=rules.canu.output.assembly,
    output:
        assembly=outdir / "{sample}" / "canu" / "{technology}" / "assembly.canu.{technology}.fasta",
        id_fofn=outdir / "{sample}" / "canu" / "{technology}" / "{sample}.contigs.nobubbles.{technology}.fofn",
    threads: 1
    resources:
        mem_mb=500
    singularity: containers["conda"]
    conda: envs["remove_bubbles"]
    params:
        pattern="'^>(?P<id>\w+)\s.*suggestBubble=no.*$'",
        replace_with="'$id'",
        extras="-uuu --no-line-number",  # disable smart filtering with -uuu
    shell:
        """
        rg {params.extras} \
            --only-matching {params.pattern} \
            --replace {params.replace_with} \
            {input.assembly} > {output.id_fofn}
        fastaq filter --ids_file {output.id_fofn} {input.assembly} {output.assembly}
        """
