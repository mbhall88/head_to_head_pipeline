"""Spades will be run on all three techs at the same time."""

rule spades:
    input:
         illumina1 = illumina_dir / "{sample}" / "{sample}.R1.fastq.gz",
         illumina2 = illumina_dir / "{sample}" / "{sample}.R2.fastq.gz",
         pacbio    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq",
         nanopore  = ont_dir / "{sample}" / "{sample}.nanopore.fastq.gz"
    output:
         assembly = outdir / "{sample}" / "spades" / "scaffolds.fasta"
    threads: 16
    resources:
             mem_mb=lambda wildcards, attempt: 32000 * attempt,
             mem_gb=lambda wildcards, attempt: 32 * attempt
    singularity: "docker://quay.io/biocontainers/spades:3.14.0--h2d02072_0"
    shell:
         """
         outdir=$(dirname {output.assembly})
         spades.py \
            -o $outdir \
            -1 {input.illumina1} \
            -2 {input.illumina2} \
            -s {input.pacbio} \
            --nanopore {input.nanopore} \
            --threads {threads} \
            --memory {resources.mem_gb} \
            --careful
         """

