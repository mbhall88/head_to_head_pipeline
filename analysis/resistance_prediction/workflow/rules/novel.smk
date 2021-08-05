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
