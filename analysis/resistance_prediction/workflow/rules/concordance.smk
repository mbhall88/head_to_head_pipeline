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


rule mutation_concordance:
    input:
        illumina_predictions=(
            RESULTS / "mykrobe/predict/illumina/{site}/{sample}/{sample}.mykrobe.json"
        ),
        nanopore_predictions=(
            RESULTS / "mykrobe/predict/nanopore/{site}/{sample}/{sample}.mykrobe.json"
        ),
    output:
        table=RESULTS / "mutation_concordance/mykrobe/{site}/{sample}.csv",
    log:
        LOGS / "mutation_concordance/{site}/{sample}.log",
    container:
        CONTAINERS["conda"]
    params:
        delim=",",
    script:
        str(SCRIPTS / "mutation_concordance.py")


rule analyse_mutation_concordance:
    input:
        tables=expand(
            RESULTS / "mutation_concordance/mykrobe/{site}/{sample}.csv",
            zip,
            site=samplesheet["site"],
            sample=samplesheet["sample"],
        ),
    output:
        RESULTS / "mutation_concordance/results.txt",
    log:
        notebook=str(RESULTS / "mutation_concordance/results.ipynb"),
    resources:
        mem_mb=int(4 * GB),
    params:
        treat_minor_as="REF",
        treat_null_as="FILT",
    conda:
        str(ENVS / "mutation_concordance.yaml")
    notebook:
        str(NOTEBOOKS / "mutation_concordance.py.ipynb")


rule analyse_results:
    input:
        coverage=QC("report/coverage.csv"),
        phenosheet=config["phenosheet"],
        concordance=all_concordance_files,
    output:
        dst_data=RESULTS / "figures/available_dst.png",
        pheno_concordance_plot=RESULTS / "figures/phenotype_concordance.png",
        geno_concordance_plot=RESULTS / "figures/illumina_concordance.png",
        pheno_concordance_csv=RESULTS / "figures/phenotype_concordance.csv",
        illumina_concordance_csv=RESULTS / "figures/illumina_concordance.csv",
        pheno_coverage_plot=RESULTS / "figures/phenotype_coverage.png",
    params:
        ignore_drugs={"pyrazinamide", "moxifloxacin"},
        unknown_is_resistant=False,
        minor_is_susceptible=True,
        failed_is_resistant=False,
    conda:
        str(ENVS / "analyse_results.yaml")
    log:
        notebook=RESULTS / "figures/analysis.processed.ipynb",
    notebook:
        str(NOTEBOOKS / "analysis.py.ipynb")
