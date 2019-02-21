rule generate_krona_input:
    input:
        bam = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam",
        bam_index = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam.bai",
        metadata = rules.download_decontamination_db.output.metadata
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


rule qc_plot:
    input:
        fastq = "analysis/{region}/nanopore/{run}/trimmed/{sample}.trimmed.fastq.gz",
        bam = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam",
        bam_index = "analysis/{region}/nanopore/{run}/mapped/{sample}.sorted.bam.bai",
    output:
        "analysis/{region}/nanopore/{run}/plotting/qc/{sample}.qc.pdf"
    threads:
        config["qc_plot"]["threads"]
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * config["qc_plot"]["memory"]
            )
    params:
        downsample = config["qc_plot"]["downsample"]
    singularity:
        config["qc_plot"]["container"]
    log:
        "analysis/logs/qc_plot_{region}_{run}_{sample}.log"
    shell:
        """
        pistis --fastq {input.fastq} \
            --output {output} \
            --bam {input.bam} \
            --downsample {params.downsample} 2> {log}
        """
