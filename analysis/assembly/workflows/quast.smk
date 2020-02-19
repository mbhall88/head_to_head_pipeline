rule quast:
    input:
        flye_pb = outdir / "{sample}" / "flye" / "pacbio" / "assembly.circularise.fasta",
        flye_pb_uc = outdir / "{sample}" / "flye" / "pacbio" / "unicycler" / "assembly.fasta",
        flye_pb_1racon = outdir / "{sample}" / "flye" / "pacbio" / "racon" / "assembly.1x.racon.fasta",
        flye_pb_1racon_pilon = outdir / "{sample}" / "flye" / "pacbio" / "racon" / "pilon" / "final.pilon.fasta",
        flye_ont = outdir / "{sample}" / "flye" / "nanopore" / "assembly.circularise.fasta",
        flye_ont_uc = outdir / "{sample}" / "flye" / "nanopore" / "unicycler" / "assembly.fasta",
        flye_ont_1racon = outdir / "{sample}" / "flye" / "nanopore" / "racon" / "assembly.1x.racon.fasta",
        flye_ont_1racon_pilon = outdir / "{sample}" / "flye" / "nanopore" / "racon" / "pilon" / "final.pilon.fasta",
        uc_pb = outdir / "{sample}" / "unicycler" / "pacbio" / "assembly.fasta",
        uc_ont = outdir / "{sample}" / "unicycler" / "nanopore" / "assembly.fasta",
        spades = outdir / "{sample}" / "spades" / "scaffolds.circularise.fasta",
        spades_pilon = outdir / "{sample}" / "spades" / "pilon" / "final.pilon.fasta",
        canu_pb = outdir / "{sample}" / "canu" / "pacbio" / "{sample}.contigs.fasta",
        canu_pb_1racon = outdir / "{sample}" / "canu" / "pacbio" / "racon" / "assembly.1x.racon.fasta",
        canu_pb_1racon_pilon = outdir / "{sample}" / "canu" / "pacbio" / "racon" / "pilon" / "final.pilon.fasta",
        canu_ont = outdir / "{sample}" / "canu" / "nanopore" / "{sample}.contigs.fasta",
        canu_ont_1racon = outdir / "{sample}" / "canu" / "nanopore" / "racon" / "assembly.1x.racon.fasta",
        canu_ont_1racon_pilon = outdir / "{sample}" / "canu" / "nanopore" / "racon" / "pilon" / "final.pilon.fasta",
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
        mem_mb = lambda wildcards, attempt: 12000 * attempt
    singularity: containers["quast"]
    params:
        labels = (
            "spades,spades_p,flye_pb,flye_pb_uc,flye_pb_1r,flye_pb_1r_p,flye_ont,"
            "flye_ont_uc,flye_ont_1r,flye_ont_1r_p,uc_pb,uc_ont,canu_pb,canu_pb_1r,"
            "canu_pb_1r_p,canu_ont,canu_ont_1r,canu_ont_1r_p"
        ),
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
             {input.spades} {input.spades_pilon} {input.flye_pb} {input.flye_pb_uc} \
             {input.flye_pb_1racon} {input.flye_pb_1racon_pilon} {input.flye_ont} \
             {input.flye_ont_uc} {input.flye_ont_1racon} {input.flye_ont_1racon_pilon} \
             {input.uc_pb} {input.uc_ont} {input.canu_pb} {input.canu_pb_1racon} \
             {input.canu_pb_1racon_pilon} {input.canu_ont} {input.canu_ont_1racon} \
             {input.canu_ont_1racon_pilon}
         """


# todo: remove some intermediate result assemblies
