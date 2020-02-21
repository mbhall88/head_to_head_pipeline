from pathlib import Path


rule quast:
    input:
        flye_pb              = outdir / "{sample}" / "flye" / "pacbio" / "assembly.flye.pacbio.fasta",
        flye_pb_polish       = outdir / "{sample}" / "flye" / "pacbio" / "polished_assembly.flye.pacbio.fasta",
        flye_ont             = outdir / "{sample}" / "flye" / "nanopore" / "assembly.flye.nanopore.fasta",
        flye_ont_polish      = outdir / "{sample}" / "flye" / "nanopore" / "polished_assembly.flye.nanopore.fasta",
        unicycler_pb         = outdir / "{sample}" / "unicycler" / "pacbio" / "assembly.unicycler.pacbio.fasta",
        unicycler_pb_polish  = outdir / "{sample}" / "unicycler" / "pacbio" / "polished_assembly.unicycler.pacbio.fasta",
        unicycler_ont        = outdir / "{sample}" / "unicycler" / "nanopore" / "assembly.unicycler.nanopore.fasta",
        unicycler_ont_polish = outdir / "{sample}" / "unicycler" / "nanopore" / "polished_assembly.unicycler.nanopore.fasta",
        spades               = outdir / "{sample}" / "spades" / "pacbio" / "assembly.spades.pacbio.fasta",
        spades_polish        = outdir / "{sample}" / "spades" / "pacbio" / "polished_assembly.spades.pacbio.fasta",
        canu_pb              = outdir / "{sample}" / "canu" / "pacbio" / "assembly.canu.pacbio.fasta",
        canu_pb_polish       = outdir / "{sample}" / "canu" / "pacbio" / "polished_assembly.canu.pacbio.fasta",
        canu_ont             = outdir / "{sample}" / "canu" / "nanopore" / "assembly.canu.nanopore.fasta",
        canu_ont_polish      = outdir / "{sample}" / "canu" / "nanopore" / "polished_assembly.canu.nanopore.fasta",
        illumina1            = rules.trim_illumina.output.forward_paired,
        illumina2            = rules.trim_illumina.output.reverse_paired,
        pacbio               = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
        nanopore             = ont_dir / "{sample}" / "{sample}.nanopore.fastq.gz",
        reference_genome     = H37RV["genome"],
        reference_features   = H37RV["features"],
    output:
        report = outdir / "{sample}" / "quast" / "report.pdf"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 12000 * attempt
    singularity: containers["quast"]
    params:
        labels = (
            "flye_pb,flye_pb_pol,flye_ont,flye_ont_pol,uc_pb,uc_pb_pol,uc_ont,"
            "uc_ont_pol,spades,spades_pol,canu_pb,canu_pb_pol,canu_ont,canu_ont_pol"
        ),
        outdir = lambda wildcards, output: Path(output.report).parent,
        extras = ""
    shell:
         """
         quast.py -o {params.outdir} \
             --threads {threads} \
             --labels  {params.labels} \
             -r {input.reference_genome} \
             --features {input.reference_features} \
             {params.extras} \
             --pe1 {input.illumina1} \
             --pe2 {input.illumina2} \
             --pacbio {input.pacbio} \
             --nanopore {input.nanopore} \
             {input.flye_pb} {input.flye_pb_polish} \
             {input.flye_ont} {input.flye_ont_polish} \
             {input.unicycler_pb} {input.unicycler_pb_polish} \
             {input.unicycler_ont} {input.unicycler_ont_polish} \
             {input.spades} {input.spades_polish} \
             {input.canu_pb} {input.canu_pb_polish} \
             {input.canu_ont} {input.canu_ont_polish}
             
         """


