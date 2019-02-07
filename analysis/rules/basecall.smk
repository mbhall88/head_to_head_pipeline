def determine_config_model(wildcards, input, output, threads, resources):
    """dna_r9.4.1_450bps.cfg"""
    pass

rule basecall:
    input:
        "data/{region}/nanopore/{run}/f5s",
    output:
        directory("analysis/{region}/nanopore/{run}/basecalling")
    threads:
        config["basecall"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["basecall"]["memory"]
    singularity:
        config["basecall"]["container"]
    params:
        model_config = determine_config_model
    log:
        "analysis/logs/basecall_{region}_{run}.log"
    shell:
        """
        guppy_basecaller --input_path {input} \
            --save_path {output} \
            --recursive \
            --worker_threads {threads} \
            --config {params.model_config} 2> {log}
        """

rule combine_fastq:
    input:
        rules.basecall.output
    output:
        temp("analysis/{region}/nanopore/{run}/basecalling/combined.fastq")
    params:
        fastqs = lambda wildcards, input: input[0] + "/*.fastq"
    log:
        "analysis/logs/combine_fastq_{region}_{run}.log"
    shell:
        """
        cat {params.fastqs} > {output} 2> {log}
        """
