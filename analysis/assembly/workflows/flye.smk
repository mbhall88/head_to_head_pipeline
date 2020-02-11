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
    singularity: config["containers"]["flye"]
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
        long_read_assembly = rules.flye.output.assembly,
    output:
        assembly = outdir / "{sample}" / "flye" / "{technology}" / "unicycler" / "assembly.fasta",
        assembly_graph = outdir / "{sample}" / "flye" / "{technology}" / "unicycler" / "assembly.gfa",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    singularity: config["containers"]["unicycler"]
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
