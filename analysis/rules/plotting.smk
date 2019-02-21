rule generate_krona_input:
    input:
        bam = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam",
        bam_index = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam.bai",
        metadata = metadata = rules.download_decontamination_db.output.metadata
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
    conda:
        "../envs/generate_krona_input.yaml"
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
