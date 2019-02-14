rule trim:
    input:
        "analysis/{region}/nanopore/{run}/demultiplex/FIX_NAMES_COMPLETE"
    output:
        "analysis/{region}/nanopore/{run}/trimmed/{sample}.trimmed.fastq.gz"
    threads:
        config["trim"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["trim"]["memory"]
    singularity:
        config["trim"]["container"]
    params:
        fastq = "analysis/{region}/nanopore/{run}/demultiplex/{sample}.fastq.gz",
    log:
        "analysis/logs/trim_{region}_{run}_{sample}.log"
    shell:
        """
        porechop --input {params.fastq}  \
            --output {output} \
            --threads {threads} \
            --format auto \
            --discard_middle > {log}
        """
