rule extract_panel_genes_from_vcf:
    input:
        annotation=RESOURCES / "h37rv.gff3",
        vcf=BuildPrg("vcfs/filtered/sparse.filtered.vcf.gz"),
        index=BuildPrg("vcfs/filtered/sparse.filtered.vcf.gz.csi"),
        oanel=RESULTS / "drprg/panel/panel.tsv",
    output:
        vcf=RESULTS / "drprg/popn_prg/popn.bcf",
    log:
        LOGS / "extract_panel_genes_from_vcf.log",
    params:
        padding=config.get("padding", 100),
        max_indel=config.get("filters", {}).get("max_indel", 20),
    conda:
        str(ENVS / "extract_panel_genes_from_vcf.yaml")
    script:
        str(SCRIPTS / "extract_panel_genes_from_vcf.py")
