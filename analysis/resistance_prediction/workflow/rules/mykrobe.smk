rule mykrobe:
    input:
        reads=QC("subsampled/{site}/{tech}/{sample}/{sample}.subsampled.fastq.gz"),
    output:
        report=RESULTS / "mykrobe/{tech}/{site}/{sample}.mykrobe.json",
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
        tech_flag=lambda wildcards: "--ont" if wildcards.tech == "nanopore" else "",
    shell:
        """
        mykrobe predict {params.flags} {params.tech_flag} --output {output.report} \
          {wildcards.sample} {params.species} > {log} 2>&1
        """
