rule concordance:
    input:
        true_pred=RESULTS / "mykrobe/illumina/{site}/{sample}.{tool}.json",
        test_pred=RESULTS / "{tool}/{tech}/{site}/{sample}.{tool}.json",
    output:
        RESULTS / "concordance/{tool}/{tech}/{site}/{sample}.{tool}.csv",
    log:
        LOGS / "concordance/{tool}/{tech}/{site}/{sample}.log",
    container:
        CONTAINERS["conda"]
    conda:
        str(ENVS / "concordance.yaml")
    params:
        opts="-r",
        script=SCRIPTS / "concordance.py",
    shell:
        """
        python {params.script} {params.opts} -a {input.true_pred} -b {input.test_pred} \
          -o {output} 2> {log}
        """
