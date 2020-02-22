from pathlib import Path


rule quast:
    input:
        flye_pb              = lambda wildcards: generate_assembly_filepath(wildcards.sample, "flye", "pacbio"),
        flye_pb_polish       = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "flye", "pacbio"),
        flye_ont             = lambda wildcards: generate_assembly_filepath(wildcards.sample, "flye", "nanopore"),
        flye_ont_polish      = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "flye", "nanopore"),
        unicycler_pb         = lambda wildcards: generate_assembly_filepath(wildcards.sample, "unicycler", "pacbio"),
        unicycler_pb_polish  = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "unicycler", "pacbio"),
        unicycler_ont        = lambda wildcards: generate_assembly_filepath(wildcards.sample, "unicycler", "nanopore"),
        unicycler_ont_polish = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "unicycler", "nanopore"),
        spades               = lambda wildcards: generate_assembly_filepath(wildcards.sample, "spades", "pacbio"),
        spades_polish        = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "spades", "pacbio"),
        canu_pb              = lambda wildcards: generate_assembly_filepath(wildcards.sample, "canu", "pacbio"),
        canu_pb_polish       = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "canu", "pacbio"),
        canu_ont             = lambda wildcards: generate_assembly_filepath(wildcards.sample, "canu", "nanopore"),
        canu_ont_polish      = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "canu", "nanopore"),
        haslr_pb             = lambda wildcards: generate_assembly_filepath(wildcards.sample, "haslr", "pacbio"),
        haslr_pb_polish      = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "haslr", "pacbio"),
        haslr_ont            = lambda wildcards: generate_assembly_filepath(wildcards.sample, "haslr", "nanopore"),
        haslr_ont_polish     = lambda wildcards: generate_polished_assembly_filepath(wildcards.sample, "haslr", "nanopore"),
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
            "haslr_pb,haslr_pb_pol,haslr_ont,haslr_ont_pol"
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
             {input.canu_ont} {input.canu_ont_polish} \
             {input.haslr_pb} {input.haslr_pb_polish} \
             {input.haslr_ont} {input.haslr_ont_polish}
         """


