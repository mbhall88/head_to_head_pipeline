rule mykrobe:
    input:
        reads=lambda wildcards: infer_reads(wildcards),
    output:
        report=(
            RESULTS / "mykrobe/predict/{tech}/{site}/{sample}/{sample}.mykrobe.json"
        ),
    shadow:
        "shallow"
    resources:
        mem_mb=int(4 * GB),
    container:
        CONTAINERS["mykrobe"]
    log:
        LOGS / "mykrobe/{tech}/{site}/{sample}.log",
    params:
        species="tb",
        flags="--force -A --debug --format json",
        tech_flag=(
            lambda wildcards: "--expected_error_rate 0.08"
            if wildcards.tech == "nanopore"
            else ""
        ),
    threads: 4
    shell:
        """
        mykrobe predict {params.tech_flag} -o {output.report} -i {input.reads} \
          {params.flags} --sample {wildcards.sample} --species {params.species} \
          -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        """
