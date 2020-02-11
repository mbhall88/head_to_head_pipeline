"""Spades will be run on all three techs at the same time."""

rule spades:
    input:
         illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
         illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
         pacbio    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
         nanopore  = ont_dir / "{sample}" / "{sample}.nanopore.fastq.gz"
    output:
         assembly = outdir / "{sample}" / "spades" / "scaffolds.fasta"
    threads: 32
    resources:
             mem_mb=lambda wildcards, attempt: 64000 + (16000 * (attempt - 1)),
             mem_gb=lambda wildcards, attempt: 64 + (16 * (attempt - 1))
    singularity: config["containers"]["spades"]
    shell:
         """
         outdir=$(dirname {output.assembly})
         spades.py \
            -o $outdir \
            --pe1-1 {input.illumina1} \
            --pe1-2 {input.illumina2} \
            --s1 {input.pacbio} \
            --nanopore {input.nanopore} \
            --threads {threads} \
            --memory {resources.mem_gb} \
            --isolate
         """

