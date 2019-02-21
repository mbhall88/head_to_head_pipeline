rule stats_pre_filter:
    input:
        "analysis/{region}/nanopore/{run}/demultiplex/FIX_NAMES_COMPLETE"
    output:
        "analysis/{region}/nanopore/{run}/stats/{sample}.pre_filter.stats.txt",
    threads:
        config["stats"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["stats"]["memory"]
            )
    params:
        fastq = "analysis/{region}/nanopore/{run}/demultiplex/{sample}.fastq.gz"
    singularity:
        config["stats"]["container"]
    log:
        "analysis/logs/stats_pre_filter_{region}_{run}_{sample}.log"
    shell:
        """
        NanoStat --fastq {params.fastq} \
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


rule report:
    input:
        qc_plot = "analysis/{region}/nanopore/{run}/plotting/qc/{sample}.qc.pdf",
        stats_pre_filter = "analysis/{region}/nanopore/{run}/stats/{sample}.pre_filter.stats.txt",
        stats_post_filter = "analysis/{region}/nanopore/{run}/stats/{sample}.post_filter.stats.txt",
        mykrobe = "analysis/{{region}}/nanopore/{{run}}/mykrobe/{{sample}}.mykrobe.{ext}".format(ext=config["mykrobe"]["output_format"]),
        trim_log = "analysis/logs/trim_{region}_{run}_{sample}.log",
    output:
        "analysis/{region}/nanopore/{run}/report/{sample}.report.html"
    threads:
        config["report"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["report"]["memory"]
            )
    params:
        sample="{sample}"
    log:
        "analysis/logs/report_{region}_{run}_{sample}.log"
    script:
        "../scripts/report.py"
