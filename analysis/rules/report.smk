rule stats_pre_filter:
    input:
        bam = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam",
        bam_index = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam.bai",
    output:
        "analysis/{region}/nanopore/{run}/stats/{sample}.pre_filter.stats.txt",
    threads:
        config["stats"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["stats"]["memory"]
            )
    singularity:
        config["stats"]["container"]
    log:
        "analysis/logs/stats_pre_filter_{region}_{run}_{sample}.log"
    shell:
        """
        NanoStat --bam {input.bam} \
            --name {output} \
            --threads {threads} 2> {log}
        """

rule stats_post_filter:
    input:
        fastq = "analysis/{region}/nanopore/{run}/trimmed/{sample}.trimmed.fastq.gz",
    output:
        "analysis/{region}/nanopore/{run}/stats/{sample}.post_filter.stats.txt"
    threads:
        config["stats"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["stats"]["memory"]
            )
    singularity:
        config["stats"]["container"]
    log:
        "analysis/logs/stats_post_filter_{region}_{run}_{sample}.log"
    shell:
        """
        NanoStat --fastq {input.fastq} \
            --name {output} \
            --threads {threads} 2> {log}
        """
