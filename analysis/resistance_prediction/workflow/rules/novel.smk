rule extract_panel_genes_from_compass_vcf:
    input:
        annotation=RESOURCES / "h37rv.gff3",
        vcf=COMPASS_DIR / "{sample}.compass.vcf.gz",
        panel=RESULTS / "drprg/panel/panel.tsv",
    output:
        vcf=RESULTS / "novel/vcfs/{sample}.res_genes.bcf",
    log:
        LOGS / "extract_panel_genes_from_compass_vcf/{sample}.log",
    params:
        padding=config.get("padding", 100),
        apply_filters=True,
        only_alt=True,
    conda:
        str(ENVS / "extract_panel_genes_from_vcf.yaml")
    script:
        str(SCRIPTS / "extract_panel_genes_from_vcf.py")


rule bcftools_index:
    input:
        rules.extract_panel_genes_from_compass_vcf.output.vcf,
    output:
        RESULTS / "novel/vcfs/{sample}.res_genes.bcf.csi",
    params:
        extra="",
    wrapper:
        "0.77.0/bio/bcftools/index"


rule subtract_panel_variants_from_compass:
    input:
        panel=rules.drprg_build.output.vcf,
        query_idx=rules.bcftools_index.output[0],
        query=rules.extract_panel_genes_from_compass_vcf.output.vcf,
    output:
        vcf=RESULTS / "novel/filtered_vcfs/{sample}.novel.vcf",
    resources:
        mem_mb=int(0.5 * GB),
    log:
        LOGS / "subtract_panel_variants_from_compass/{sample}.log",
    conda:
        str(ENVS / "subtract_variants.yaml")
    script:
        str(SCRIPTS / "subtract_variants.py")


rule assess_drprg_novel_calls:
    input:
        truth_asm=rules.drprg_build.output.ref,
        vcf_ref=rules.drprg_build.output.ref,
        vcf_to_eval=rules.subtract_panel_variants_from_drprg.output.vcf,
        truth_vcf=rules.subtract_panel_variants_from_compass.output.vcf,
    output:
        summary=RESULTS / "novel/assessment/{tech}/{site}/{sample}/summary_stats.json",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: int(4 * GB) * attempt,
    params:
        options="--force --filter_pass PASS,.",
        flank_length=100,
        outdir=lambda wildcards, output: Path(output.summary).parent,
    log:
        LOGS / "assess_drprg_novel_calls/{tech}/{site}/{sample}.log",
    conda:
        str(ENVS / "varifier.yaml")
    shell:
        """
        rm -rf {params.outdir}
        varifier vcf_eval {params.options} \
            --flank_length {params.flank_length} \
            --truth_vcf {input.truth_vcf} \
            {input.truth_asm} \
            {input.vcf_ref} \
            {input.vcf_to_eval} \
            {params.outdir} 2> {log}
        """
