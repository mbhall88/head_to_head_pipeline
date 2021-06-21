rule mykrobe:
    input:
        reads=lambda wildcards: infer_reads(wildcards),
    output:
        report=RESULTS / "mykrobe/predict/{tech}/{site}/{sample}/{sample}.mykrobe.json",
    shadow:
        "shallow"
    resources:
        mem_mb=int(2 * GB),
    container:
        CONTAINERS["mykrobe"]
    log:
        LOGS / "mykrobe/{tech}/{site}/{sample}.log",
    params:
        species="tb",
        flags="--force --report_all_calls --debug --format json",
        tech_flag=(
            lambda wildcards: "--expected_error_rate 0.08 --ploidy haploid"
            if wildcards.tech == "nanopore"
            else ""
        ),
    shell:
        """
        mykrobe predict {params.tech_flag} --output {output.report} --seq {input.reads} \
          {params.flags} {wildcards.sample} {params.species} > {log} 2>&1
        """
