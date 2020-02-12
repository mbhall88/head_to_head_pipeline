rule quast:
    input:
        flye_pb = outdir / "{sample}" / "flye" / "pacbio" / "assembly.fasta",
        flye_pb_uc = outdir / "{sample}" / "flye" / "pacbio" / "unicycler" / "assembly.fasta",
        flye_pb_1racon = outdir / "{sample}" / "flye" / "pacbio" / "racon" / "assembly.1x.racon.fasta",
        flye_ont = outdir / "{sample}" / "flye" / "nanopore" / "assembly.fasta",
        flye_ont_uc = outdir / "{sample}" / "flye" / "nanopore" / "unicycler" / "assembly.fasta",
        flye_ont_1racon = outdir / "{sample}" / "flye" / "nanopore" / "racon" / "assembly.1x.racon.fasta",
        spades = outdir / "{sample}" / "spades" / "scaffolds.fasta",
        illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
        illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
        pacbio    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
        nanopore  = ont_dir / "{sample}" / "{sample}.nanopore.fastq.gz",
        reference_genome = H37RV["genome"],
        reference_features = H37RV["features"],
    output:
        report = outdir / "{sample}" / "quast" / "report.pdf"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt
    singularity: containers["quast"]
    params:
        labels = "spades,flye_pb,flye_pb_uc,flye_pb_1r,flye_ont,flye_ont_uc,flye_ont_1r",
        extras = ""
    shell:
         """
         outdir=$(dirname {output.report})
         quast.py -o $outdir \
             --threads {threads} \
             --labels  {params.labels} \
             -r {input.reference_genome} \
             --features {input.reference_features} \
             {params.extras} \
             --pe1 {input.illumina1} \
             --pe2 {input.illumina2} \
             --pacbio {input.pacbio} \
             --nanopore {input.nanopore} \
             {input.spades} {input.flye_pb} {input.flye_pb_uc} {input.flye_pb_1racon} \
             {input.flye_ont} {input.flye_ont_uc} {input.flye_ont_1racon}
         """
