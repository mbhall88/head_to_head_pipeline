rule quast:
    input:
        flye_pb = outdir / "{sample}" / "flye" / "pacbio" / "assembly.fasta",
        flye_ont = outdir / "{sample}" / "flye" / "nanopore" / "assembly.fasta",
        spades = outdir / "{sample}" / "spades" / "scaffolds.fasta",
        illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
        illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
        pacbio    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
        nanopore  = ont_dir / "{sample}" / "{sample}.nanopore.fastq.gz",
    output:
        report = outdir / "{sample}" / "quast" / "report.pdf"
    threads = 8
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    params:
        genome_size = config["genome_size"]
    shell:
        """
        outdir=$(dirname {output.report})
        quast.py -o $outdir \
            --threads {threads} \
            --labels spades,flye_pb,flye_ont \
            --gene-finding \
            --conserved-genes-finding \
            --est-ref-size {params.genome_size} \
            --pe1 {input.illumina1} \
            --pe2 {input.illumina2} \
            --pacbio {input.pacbio} \
            --nanopore {input.nanopore} \
            {input.spades} {input.flye_pb} {input.flye_ont}
        """
