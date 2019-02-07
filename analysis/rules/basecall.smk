 class InvalidFlowcell(Exception):
     __module__ = Exception.__module__


def determine_config_model(wildcards, input, output, threads, resources):
    """dna_r9.4.1_450bps.cfg"""
    run_flowcell = set(samples.xs((wildcards.region, wildcards.run)).loc[ : , 'flowcell_version'])
    
    if len(run_flowcell) != 1:
        raise InvalidFlowcell("Multiple flowcell versions found for the same run - {}".format(wildcards.run))

    flowcell_ver = next(run_flowcell)

    if flowcell_ver.starswith("R9.5"):
        return "dna_r9.5_450bps.cfg"
    elif flowcell.startswith("R9.4"):
        return "dna_r9.4.1_450bps.cfg" 
    else:
        raise InvalidFlowcell("Currently only support flowcells beginning with R9.4 or R9.5. You provided {}".format(flowcell_ver))
    

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
