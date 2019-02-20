rule generate_krona_input:
    input:
        bam = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam",
        bam_index = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam.bai",
    output:
        temp("analysis/{region}/nanopore/{run}/plotting/{sample}.krona.txt")
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
