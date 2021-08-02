rule concordance:
    input:
        true_pred=(
            RESULTS / "mykrobe/predict/illumina/{site}/{sample}/{sample}.mykrobe.json"
        ),
        test_pred=(
            RESULTS / "{tool}/predict/{tech}/{site}/{sample}/{sample}.{tool}.json"
        ),
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


rule analyse_results:
    input:
        coverage=QC("report/coverage.csv"),
        phenosheet=config["phenosheet"],
        concordance=[
            RESULTS / f"concordance/{tool}/{tech}/{site}/{sample}.{tool}.csv"
            for site, sample, tool, tech in zip(
                samplesheet["site"], samplesheet["sample"], TOOLS, TECHS
            )
        ],
    output:
        dst_data=RESULTS / "figures/available_dst.png",
        pheno_concordance_plot=RESULTS / "figures/phenotype_concordance.png",
        pheno_concordance_csv=RESULTS / "figures/phenotype_concordance.csv",
        illumina_concordance_csv=RESULTS / "figures/illumina_concordance.csv",
        pheno_coverage_plot=RESULTS / "figures/phenotype_coverage.png",
    params:
        ignore_drugs={"pyrazinamide", "moxifloxacin"},
        unknown_is_resistant=False,
        minor_is_susceptible=False,
    conda:
        str(ENVS / "analyse_results.yaml")
    log:
        notebook=RESULTS / "figures/analysis.processed.ipynb",
    notebook:
        str(NOTEBOOKS / "analysis.py.ipynb")
