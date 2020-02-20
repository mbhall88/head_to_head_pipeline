"""Spades will be run on all three techs at the same time."""
from pathlib import Path

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
    singularity: containers["spades"]
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

"""We remove '.circularised' from the output headers as they cause issues in prokka."""
rule circularise_spades:
    input:
         assembly = rules.spades.output.assembly,
    output:
          assembly = outdir / "{sample}" / "spades" / "scaffolds.circularise.fasta",
    threads: 1
    resources:
             mem_mb = lambda wildcards, attempt: 1000 * attempt
    params:
          output_prefix = lambda wildcards, input: Path(input.assembly).with_suffix("")
    singularity: containers["circlator"]
    shell:
         """
         circlator minimus2 {input.assembly} {params.output_prefix}
         sed -ir 's/\.circularised//g' {output.assembly}
         """


rule pilon_polish_spades:
    input:
         assembly = rules.circularise_spades.output.assembly,
         illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
         illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
    output:
          polished_assembly = outdir / "{sample}" / "spades" / "pilon" / "final.pilon.fasta"
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


rule annotate_spades:
    input:
        assembly = rules.pilon_polish_spades.output.polished_assembly,
    output:
        annotation = outdir / "{sample}" / "spades" / "prokka" / "spades.pilon.gff",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    singularity: containers["prokka"]
    params:
        annotation = config["h37rv"]["annotation"],
        outdir = lambda wildcards, output: Path(output.annotation).parent,
        prefix = lambda wildcards, output: Path(output.annotation).with_suffix("").name,
        extras = "--force"
    shell:
        """
        prokka --proteins {params.annotation} \
            --cpus {threads} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            {params.extras} \
            {input.assembly}
        """


rule map_illumina_reads_to_spades_polished_assembly:
    input:
         assembly = rules.pilon_polish_spades.output.polished_assembly,
         illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
         illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
    output:
          bam = outdir / "{sample}" / "spades" / "mapping" / "{sample}.pilon.illumina.spades.bam",
          bam_index = outdir / "{sample}" / "spades" / "mapping" / "{sample}.pilon.illumina.spades.bam.bai",
    threads: 8
    resources:
             mem_mb = lambda wildcards, attempt: 8000 * attempt
    singularity: containers["conda"]
    conda: envs["aln_tools"]
    params:
        view_extras = "-hu", # include header, uncompressed BAM
        filter_flags = "-F {}".format(config["illumina_flag_filter"]),
    shell:
         """
         bwa index {input.assembly}
         bwa mem -t {threads} {input.assembly} {input.illumina1} {input.illumina2} | \
             samtools view {params.filter_flags} {params.view_extras} - | \
             samtools sort -@ {threads} -o {output.bam} -
         samtools index {output.bam}
         """

rule stats_spades:
    input:
        bam = rules.map_illumina_reads_to_spades_polished_assembly.output.bam,
    output:
        stats = outdir / "{sample}" / "spades" / "qc" / "{sample}.pilon.illumina.stats"
    threads: 1
    resources:
        mem_mb = 1000
    singularity: containers["samtools"]
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        """

rule mpileup_spades:
    input:
        bam = rules.map_illumina_reads_to_spades_polished_assembly.output.bam,
        assembly = rules.pilon_polish_spades.output.polished_assembly,
    output:
        pileup = outdir / "{sample}" / "spades" / "mapping" / "{sample}.pilon.illumina.spades.pileup",
    singularity: containers["samtools"]
    params:
        extras = "-aa"
    shell:
        """
        samtools mpileup {params.extras} \
            -o {output.pileup} \
            --fasta-ref {input.assembly} \
            {input.bam}
        """

rule assess_per_base_accuracy_spades:
    input:
        bam = rules.map_illumina_reads_to_spades_polished_assembly.output.bam,
        pileup = rules.mpileup_spades.output.pileup,
    output:
        json = outdir / "{sample}" / "spades" / "assessment" / "{sample}.spades.pilon.json",
        bed = outdir / "{sample}" / "spades" / "assessment" / "{sample}.spades.pilon.bed"
    params:
        script = scripts["assess_per_base"],
        min_depth = config["min_depth"],
        quorum = config["quorum"],
        prefix = lambda wildcards, output: Path(output.json).with_suffix(""),
    singularity: containers["conda"]
    conda: envs["assess_per_base"]
    shell:
        """
        python3 {params.script} \
            --bam {input.bam} \
            --pileup {input.pileup} \
            --min-depth {params.min_depth} \
            --quorum {params.quorum} \
            --prefix {params.prefix}
        """
