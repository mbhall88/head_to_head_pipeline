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
        assembly      = outdir / "{sample}" / "flye" / "{technology}" / "assembly.fasta",
        assembly_info = outdir / "{sample}" / "flye" / "{technology}" / "assembly_info.txt",
        assembly_graph = outdir / "{sample}" / "flye" / "{technology}" / "assembly_graph.gfa",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 16000 * attempt,
    singularity: containers["flye"]
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


"""
This rule is effectively using Unicycler to polish the Flye assembly.
"""
rule unicycler_polish_flye:
    input:
        illumina1          = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
        illumina2          = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
        long_reads         = mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
        long_read_assembly = rules.flye.output.assembly_graph,
    output:
        assembly = outdir / "{sample}" / "flye" / "{technology}" / "unicycler" / "assembly.fasta",
        assembly_graph = outdir / "{sample}" / "flye" / "{technology}" / "unicycler" / "assembly.gfa",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 16000 * attempt
    singularity: containers["unicycler"]
    params:
        verbosity = 2,
    shell:
        """
        outdir=$(dirname {output.assembly})
        unicycler --short1 {input.illumina1} \
            --short2 {input.illumina2} \
            --long {input.long_reads} \
            --out $outdir \
            --verbosity {params.verbosity} \
            --threads {threads} \
            --existing_long_read_assembly {input.long_read_assembly}
        """

rule circularise_flye:
    input:
        assembly = rules.flye.output.assembly,
    output:
        assembly = outdir / "{sample}" / "flye" / "{technology}" / "assembly.circularise.fasta",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1000 * attempt
    params:
        output_prefix = lambda wildcards, input: Path(input.assembly).with_suffix("")
    singularity: containers["circlator"]
    shell:
        """
        circlator minimus2 {input.assembly} {params.output_prefix}
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


rule map_long_reads_to_flye_assembly:
    input:
        reads = mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
        flye_assembly = rules.circularise_flye.output.assembly,
    output:
        sam = outdir / "{sample}" / "flye" / "{technology}" / "mapping" / "{sample}.{technology}.flye.sam"
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
            {input.flye_assembly} \
            {input.reads} > {output.sam}
        """

rule map_illumina_reads_to_flye_assembly:
    input:
        assembly = rules.circularise_flye.output.assembly,
        illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
        illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
    output:
        bam = outdir / "{sample}" / "flye" / "{technology}" / "mapping" / "{sample}.illumina.flye.bam",
        bam_index = outdir / "{sample}" / "flye" / "{technology}" / "mapping" / "{sample}.illumina.flye.bam.bai",
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    shell:
        """
        bwa index {input.assembly}
        bwa mem -t {threads} {input.assembly} {input.illumina1} {input.illumina2} | \
            samtools sort -@ {threads} -o {output.bam} -
        samtoools index {output.bam}
        """


rule racon_polish_flye:
    input:
        reads = mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
        sam = rules.map_long_reads_to_flye_assembly.output.sam,
        assembly = rules.circularise_flye.output.assembly
    output:
        polished_assembly = outdir / "{sample}" / "flye" / "{technology}" / "racon" / "assembly.1x.racon.fasta"
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
rule pilon_polish_flye:
    input:
        assembly = rules.racon_polish_flye.output.polished_assembly,
        illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
        illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
    output:
          polished_assembly = outdir / "{sample}" / "flye" / "{technology}" / "racon" / "pilon" / "final.pilon.fasta"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    params:
        pilon_url = f"https://github.com/broadinstitute/pilon/releases/download/v{config['pilon_version']}/pilon-{config['pilon_version']}.jar",
        script_url = "https://raw.githubusercontent.com/mbhall88/bioscripts/master/python/pilon_iterative.py",
        mem_gb = lambda wildcards, resources: int(resources.mem_mb / 1000),
        max_iterations = 10,
        outdir = lambda wildcards, output: Path(output.polished_assembly).parent,
        final_fasta = lambda wildcards, output: Path(output.polished_assembly).name,
    shell:
         """
         script=$(realpath pilon.py)
         wget -O $script {params.script_url}
         pilon_jar=$(realpath pilon.jar)
         wget -O $pilon_jar {params.pilon_url}
         
         bwa index {input.assembly}
         python3 $script \
             --pilon_java_xmx {params.mem_gb}G \
             --threads {threads} \
             --max_iterations {params.max_iterations} \
             --pilon_jar $pilon_jar \
             --assembly_fasta {input.assembly} \
             --reads1 {input.illumina1} \
             --reads2 {input.illumina2} \
             --outdir {params.outdir} \
             --final_fasta {params.final_fasta}
         rm $pilon_jar $script
         """

#todo: add hypo rule - polish with CCS and ONT and then with Illumina or vice versa..
