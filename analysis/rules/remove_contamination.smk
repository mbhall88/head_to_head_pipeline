rule download_decontamination_db:
    output:
        fasta    = "data/decontamination_db/remove_contam.fa.gz",
        metadata = "data/decontamination_db/remove_contam.tsv"
    threads:
        config["download_decontamination_db"]["memory"]
    resources:
        mem_mb = (
    lambda wildcards, attempt: attempt * config["download_decontamination_db"]["memory"]
            )
    params:
        script  = config["download_decontamination_db"]["perl_script"],
        out_dir = config["download_decontamination_db"]["out_dir"]
    log:
        "analysis/logs/download_decontamination_db.log"
    shell:
        """
        bash analysis/scripts/make_decontamination_db.sh {params.script} \
            {params.out_dir} \
            {output.fasta} \
            output.metadata 2> {log}
        """

rule index_decontamination_db:
    input:
        rules.download_decontamination_db.output.fasta
    output:
        index = "data/decontamination_db/remove_contam.{}.mmi".format(config["index_decontamination_db"]["preset"])
    threads: 1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["index_decontamination_db"]["memory"]
            )
    params:
        preset = config["index_decontamination_db"]["preset"],
        split  = config["index_decontamination_db"]["split"]
    singularity:
        config["index_decontamination_db"]["container"]
    log:
        "analysis/logs/index_decontamination_db.log"
    shell:
        """
        minimap2 -I {params.split} \
            -x {params.preset} \
            -d {output.index} \
            {input} 2> {log}
        """

rule mapping:
    input:
        index = rules.index_decontamination_db.output.index,
        query = "analysis/{region}/nanopore/{run}/trimmed/{sample}.trimmed.fastq.gz"
    output:
        temp("analysis/{region}/nanopore/{run}/mapped/{sample}.sam")
    threads:
        config["mapping"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["mapping"]["memory"]
            )
    params:
        preset = config["index_decontamination_db"]["preset"],
    singularity:
        config["mapping"]["container"]
    log:
        "analysis/logs/mapping_{region}_{run}_{sample}.log"
    shell:
        """
        minimap2 -t {threads} \
            -ax {params.preset} \
            {input.index} \
            {input.query} > {output} 2> {log}
        """


rule sort:
    input:
        "analysis/{region}/nanopore/{run}/mapped/{sample}.sam"
    output:
        "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam"
    threads:
        config["sort"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["sort"]["memory"]
            )
    singularity:
        config["sort"]["container"]
    log:
        "analysis/logs/sort_{region}_{run}_{sample}.log"
    shell:
        "samtools sort -@ {threads} {input} 2> {log} > {output}"


rule index_bams:
    input:
        "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam"
    output:
        "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam.bai"
    threads:
        config["index_bams"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["index_bams"]["memory"]
            )
    singularity:
        config["index_bams"]["container"]
    log:
        "analysis/logs/index_bams_{region}_{run}_{sample}.log"
    shell:
        "samtools index -b {input} 2> {log}"


rule filter:
    input:
        bam = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam",
        index = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam.bai",
        metadata = rules.download_decontamination_db.output.metadata,
    output:
        fastq = "analysis/{region}/nanopore/{run}/filtered/{sample}.filtered.fastq.gz"
    threads:
        config["filter"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["filter"]["memory"]
            )
    singularity:
        config["filter"]["container"]
    conda:
        config["filter"]["env"]
    log:
        "analysis/logs/filter_{region}_{run}_{sample}.log"
    script:
        config["filter"]["script"]
