rule concordance:
    input:
        true_pred=RESULTS / "{tool}/illumina/{site}/{sample}.{tool}.json",
        test_pred=RESULTS / "{tool}/nanopore/{site}/{sample}.{tool}.json",
    output:
        RESULTS / "concordance/{tool}/{site}/{sample}.{tool}.csv",
    log:
        LOGS / "concordance/{tool}/{site}/{sample}.log",
    container:
        CONTAINERS["conda"]
    conda:
        ENVS / "concordance.yaml"
    params:
        opts="-r",
        script=SCRIPTS / "concordance.py",
    shell:
        """
        python {params.script} {params.opts} -a {input.true_pred} -b {input.test_pred} \
          -o {output} 2> {log}
        """
