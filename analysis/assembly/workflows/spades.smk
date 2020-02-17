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
         """


rule map_illumina_reads_to_spades_assembly:
    input:
         assembly = rules.circularise_spades.output.assembly,
         illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
         illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
    output:
          bam = outdir / "{sample}" / "spades" / "mapping" / "{sample}.illumina.spades.bam",
          bam_index = outdir / "{sample}" / "spades" / "mapping" / "{sample}.illumina.spades.bam.bai",
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



"""This rule assumes bwa mem and java are available on PATH"""
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
# todo: add rule to analyse pileup
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
    shell:
         """
         bwa index {input.assembly}
         bwa mem -t {threads} {input.assembly} {input.illumina1} {input.illumina2} | \
             samtools sort -@ {threads} -o {output.bam} -
         samtoools index {output.bam}
         """

rule assess_spades:
    input:
        assembly = rules.pilon_polish_spades.output.polished_assembly,
        bam = rules.map_illumina_reads_to_spades_polished_assembly.output.bam,
    output:
        bed = outdir / "{sample}" / "spades" / "assessment" / "spades.pilon.bed",
        positions_masked = outdir / "{sample}" / "spades" / "assessment" / "spades.pilon.txt",
    params:
        bam_to_low_qual_mask_script = config["bam_to_low_qual_mask_script"],
        depth_cutoff = config["depth_cutoff_mask"],
        assess_per_base_script = config["assess_per_base_script"],

    shell:
        """
        perl {params.bam_to_low_qual_mask_script} {params.depth_cutoff} {input.bam} {output.bed}
        python3 {params.assess_per_base_script} --bed {output.bed} --assembly {input.assembly} &> {output.positions_masked}
        sum qualities ...
        """
# todo: sum quality scores