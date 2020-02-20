from pathlib import Path
from enum import Enum

class CanuInputType(Enum):
    PACBIO = "-pacbio-hifi"
    NANOPORE = "-nanopore-raw"


def infer_canu_input_type(technology: str) -> str:
    return CanuInputType[technology.upper()].value

rule canu:
    input:
        reads = mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
        assembly = outdir / "{sample}" / "canu" / "{technology}" / "{sample}.contigs.fasta",
        assembly_graph = outdir / "{sample}" / "canu" / "{technology}" / "{sample}.unitigs.gfa",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    params:
        min_read_length = 1000,
        min_overlap_length = 500,
        genome_size = config["genome_size"],
        input_type = lambda wildcards: infer_canu_input_type(wildcards.technology),
        mem_gb = lambda wildcards, resources: int(resources.mem_mb / 1000),
        outprefix = lambda wildcards, output: Path(output.assembly).name.replace(".contigs.fasta", ""),
        outdir = lambda wildcards, output: Path(output.assembly).parent,
        extra = "corMaxEvidenceErate=0.15",
        #can try corMaxEvidenceErate=0.15 for GC rich, repetitive genomes from docs https://canu.readthedocs.io/en/latest/faq.html#my-genome-is-at-or-gc-rich-do-i-need-to-adjust-parameters-what-about-highly-repetitive-genomes
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
        """


rule remove_bubbles_canu:
    input:
        assembly = rules.canu.output.assembly,
    output:
        assembly = outdir / "{sample}" / "canu" / "{technology}" / "{sample}.contigs.nobubbles.fasta",
        id_fofn = outdir / "{sample}" / "canu" / "{technology}" / "{sample}.contigs.nobubbles.fofn",
    threads: 1
    resources:
        mem_mb = 500
    singularity: containers["fastaq"]
    params:
        pattern = "'s/^>\(\S*\)\s.*suggestBubble=no.*$/\1/p'", # thank you https://unix.stackexchange.com/a/278377/318480
        extras = "-n",
    shell:
        """
        sed {params.pattern} {input.assembly} > {output.id_fofn}
        fastaq filter --ids_file {output.id_fofn} {input.assembly} {output.assembly}
        """

"""
Taken from https://www.biorxiv.org/content/biorxiv/early/2019/08/13/635037.full.pdf
Polishing PacBio CCS with racon
minimap2 -ax map-pb --eqx -m 5000 -t {threads} --secondary=no {ref} {fastq} | samtools view -F 1796-> {sam}
racon {fastq} {sam} {ref} -u -t {threads} > {output.fasta}
"""
class MinimapPresets(Enum):
    PACBIO = "asm20"
    NANOPORE = "map-ont"


def infer_minimap2_preset(technology: str) -> str:
    return MinimapPresets[technology.upper()].value


rule map_pacbio_reads_to_canu_assembly:
    input:
        reads    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
        canu_assembly = rules.remove_bubbles_canu.output.assembly,
    output:
        sam = outdir / "{sample}" / "canu" / "{technology}" / "mapping" / "{sample}.canu.sam"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    singularity: containers["minimap2"]
    params:
        preset = lambda wildcards: infer_minimap2_preset(wildcards.technology),
        extras = "--secondary=no"
    shell:
        """
        minimap2 -ax {params.preset} \
            {params.extras} \
            -t {threads} \
            {input.canu_assembly} \
            {input.reads} > {output.sam}
        """

"""Only polish with pacbio reads"""
rule racon_polish_canu:
    input:
        reads    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
        sam = outdir / "{sample}" / "canu" / "{technology}" / "mapping" / "{sample}.canu.sam",
        assembly = rules.remove_bubbles_canu.output.assembly
    output:
        polished_assembly = outdir / "{sample}" / "canu" / "{technology}" / "racon" / "assembly.1x.racon.fasta"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    singularity: containers["racon"]
    params:
        extras = "--include-unpolished"
    shell:
        """
        racon --threads {threads} \
            {params.extras} \
            {input.reads} \
            {input.sam} \
            {input.assembly} > {output.polished_assembly}
        """

"""This rule assumes bwa mem and java are available on PATH"""
rule pilon_polish_canu:
    input:
        assembly = rules.racon_polish_canu.output.polished_assembly,
        illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
        illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
    output:
        polished_assembly = outdir / "{sample}" / "canu" / "{technology}" / "racon" / "pilon" / "final.pilon.fasta"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    params:
          pilon_jar = scripts["pilon_jar"],
          script = scripts["pilon"],
          mem_gb = lambda wildcards, resources: int(resources.mem_mb / 1000),
          max_iterations = 10,
          outdir = lambda wildcards, output: Path(output.polished_assembly).parent,
          final_fasta = lambda wildcards, output: Path(output.polished_assembly).name,
    singularity: containers["conda"]
    conda: envs["pilon"]
    shell:
         """
         bwa index {input.assembly}
         python3 {params.script} \
             --pilon_java_xmx {params.mem_gb}G \
             --threads {threads} \
             --max_iterations {params.max_iterations} \
             --pilon_jar {params.pilon_jar} \
             --assembly_fasta {input.assembly} \
             --reads1 {input.illumina1} \
             --reads2 {input.illumina2} \
             --outdir {params.outdir} \
             --final_fasta {params.final_fasta}
         """

# todo: add rule to annotate
# todo: add rule to analyse pileup
