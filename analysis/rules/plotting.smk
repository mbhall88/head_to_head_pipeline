rule create_taxonomy_lookup:
    input:
        metadata = rules.download_decontamination_db.output.metadata
    output:
        taxonomy = "data/decontamination_db/taxonomy.json"
    threads:
        config["create_taxonomy_lookup"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["create_taxonomy_lookup"]["memory"]
            )
    params:
        email = config["create_taxonomy_lookup"]["email"],
        api_key = config["create_taxonomy_lookup"]["api_key"],
    singularity:
        config["create_taxonomy_lookup"]["container"]
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/create_taxonomy_lookup.py"

rule generate_krona_input:
    input:
        bam = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam",
        bam_index = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam.bai",
        taxonomy = rules.create_taxonomy_lookup.output.taxonomy
    output:
        "analysis/{region}/nanopore/{run}/plotting/krona/{sample}.krona.txt"
    threads:
        config["generate_krona_input"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["generate_krona_input"]["memory"]
            )
    singularity:
        config["generate_krona_input"]["container"]
    script:
        "../scripts/generate_krona_input.py"

rule plot_sample_composition:
    input:
        "analysis/{region}/nanopore/{run}/plotting/krona/{sample}.krona.txt"
    output:
        "analysis/{region}/nanopore/{run}/plotting/krona/{sample}.krona.html"
    threads:
        config["plot_sample_composition"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["plot_sample_composition"]["memory"]
            )
    singularity:
        config["plot_sample_composition"]["container"]
    log:
        "analysis/logs/plot_sample_composition_{region}_{run}_{sample}.log"
    shell:
        """
        ktImportText {input} -o {output} 2> {log}
        """
