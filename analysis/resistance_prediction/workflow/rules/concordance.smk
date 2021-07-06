rule concordance:
    input:
        true_pred=RESULTS / "mykrobe/predict/illumina/{site}/{sample}/{sample}.mykrobe.json",
        test_pred=RESULTS / "{tool}/predict/{tech}/{site}/{sample}/{sample}.{tool}.json",
    output:
        RESULTS / "concordance/{tool}/{tech}/{site}/{sample}.{tool}.csv",
    log:
        LOGS / "concordance/{tool}/{tech}/{site}/{sample}.log",
    conda:
        str(ENVS / "concordance.yaml")
    params:
        opts="",
        script=SCRIPTS / "concordance.py",
    shell:
        """
        python {params.script} {params.opts} -a {input.true_pred} -b {input.test_pred} \
          -o {output} 2> {log}
        """
