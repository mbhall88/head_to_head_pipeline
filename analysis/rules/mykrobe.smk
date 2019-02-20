rule mykrobe:
    input:
        "analysis/{region}/nanopore/{run}/filtered/{sample}.filtered.fastq.gz"
    output:
        "analysis/{{region}}/nanopore/{{run}}/mykrobe/{{sample}}.mykrobe.{ext}".format(ext=config["mykrobe"]["output_format"])
    threads:
        config["mykrobe"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["mykrobe"]["memory"]
            )
    params:
        species = config["mykrobe"]["species"],
        tmp_dir = "analysis/{region}/nanopore/{run}/mykrobe/tmp_dirs/{sample}",
        output_format = config["mykrobe"]["output_format"],
    singularity:
        config["mykrobe"]["container"]
    log:
        "analysis/logs/mykrobe_{region}_{run}_{sample}.log"
    shell:
        """
        mykrobe predict {wildcards.sample} \
            {params.species} \
            --tmp {params.tmp_dir} \
            --ont \
            --seq {input} \
            --output {output} \
            --format {params.output_format}
        """
