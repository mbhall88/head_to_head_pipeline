from pathlib import Path


rule canu:
    input:
         pacbio    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
         nanopore  = ont_dir / "{sample}" / "{sample}.nanopore.fastq.gz",
    output:
        assembly = outdir / "{sample}" / "canu" / "{sample}.contigs.fasta",
        assembly_graph = outdir / "{sample}" / "canu" / "{sample}.unitigs.gfa",
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    params:
        min_read_length = 700,
        min_overlap_length = 400,
        genome_size = config["genome_size"],
        mem_gb = lambda wildcards, resources: int(resources.mem_mb / 1000),
        outprefix = lambda wildcards, output: Path(output.assembly).name.replace(".contigs.fasta", ""),
        outdir = lambda wildcards, output: Path(output.assembly).parent,
        extra = "",
        #can try corMaxEvidenceErate=0.15 for GC rich, repetitive genomes from docs https://canu.readthedocs.io/en/latest/faq.html#my-genome-is-at-or-gc-rich-do-i-need-to-adjust-parameters-what-about-highly-repetitive-genomes
    singularity: containers["canu"]
    shell:
        """
        canu -p {params.outprefix} \
            -d {params.outdir} \
            minReadLength={params.min_read_length} \
            minOverlapLength={params.min_overlap_length} \
            genomeSize={params.genome_size} \
            minThreads={threads} maxThreads={threads} \
            maxMemory={params.mem_gb} \
            -nanopore-raw {input.nanopore} \
            -pacbio-hifi {input.pacbio} {params.extra}
        """

rule circularise_canu:
    input:
         assembly = rules.canu.output.assembly,
    output:
          assembly = outdir / "{sample}" / "canu" / "{sample}.contigs.circularise.fasta",
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

"""
Taken from https://www.biorxiv.org/content/biorxiv/early/2019/08/13/635037.full.pdf
Polishing PacBio CCS with racon
minimap2 -ax map-pb --eqx -m 5000 -t {threads} --secondary=no {ref} {fastq} | samtools view -F 1796-> {sam}
racon {fastq} {sam} {ref} -u -t {threads} > {output.fasta}
"""
rule map_long_reads_to_canu_assembly:
    input:
        reads    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
        canu_assembly = rules.circularise_canu.output.assembly,
    output:
        sam = outdir / "{sample}" / "canu" / "mapping" / "{sample}.canu.sam"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 16000 * attempt
    singularity: containers["minimap2"]
    params:
        preset = "asm20",
        extras = "--secondary=no"
    shell:
        """
        minimap2 -ax {params.preset} \
            {params.extras} \
            -t {threads} \
            {input.canu_assembly} \
            {input.reads} > {output.sam}
        """

rule racon_polish_canu:
    input:
        reads    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
        sam = rules.map_long_reads_to_canu_assembly.output.sam,
        assembly = rules.circularise_canu.output.assembly
    output:
        polished_assembly = outdir / "{sample}" / "canu" / "racon" / "assembly.1x.racon.fasta"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    singularity: containers["racon"]
    params:
        extras = "--include-unpolished"
    shell:
        """
        racon --threads {threads} \
            {params.extras} \
            {input.reads} \
            {input.sam} \
            {input.assembly} > {output.polished_assembly}
        """

"""This rule assumes bwa mem and java are available on PATH"""
rule pilon_polish_canu:
    input:
        assembly = rules.racon_polish_canu.output.polished_assembly,
        illumina1 = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
        illumina2 = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
    output:
        polished_assembly = outdir / "{sample}" / "canu" / "racon" / "pilon" / "final.pilon.fasta"
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
