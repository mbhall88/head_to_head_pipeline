rule trim:
    input:
        "analysis/{region}/nanopore/{run}/demultiplex/FIX_NAMES_COMPLETE",
        fastq = "analysis/{region}/nanopore/{run}/demultiplex/{sample}.fastq.gz"
    output:
        "analysis/{region}/nanopore/{run}/trimmed/{sample}.trimmed.fastq.gz"
    threads:
        config["trim"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["trim"]["memory"]
    singularity:
        config["trim"]["container"]
    params:
        "--discard_middle ",
        "--format auto "
    log:
        "analysis/logs/trim_{region}_{run}_{sample}.log"
    shell:
        """
        porechop --input {input}  \
            --output {output} \
            --threads {threads} \
            {params} > {log}
        """
