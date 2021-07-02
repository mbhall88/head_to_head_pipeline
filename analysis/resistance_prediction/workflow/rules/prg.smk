rule extract_panel_genes_from_vcf:
    input:
        annotation=RESOURCES / "h37rv.gff3",
        vcf=BuildPrg("vcfs/filtered/sparse.filtered.vcf.gz"),
        index=BuildPrg("vcfs/filtered/sparse.filtered.vcf.gz.csi"),
        panel=RESULTS / "drprg/panel/panel.tsv",
    output:
        vcf=RESULTS / "drprg/popn_prg/popn.bcf",
    log:
        LOGS / "extract_panel_genes_from_vcf.log",
    params:
        padding=config.get("padding", 100),
    conda:
        str(ENVS / "extract_panel_genes_from_vcf.yaml")
    script:
        str(SCRIPTS / "extract_panel_genes_from_vcf.py")


rule create_references:
    input:
        genome=RESOURCES / "h37rv.fa",
        faidx=RESOURCES / "h37rv.fa.fai",
        annotation=rules.extract_panel_genes_from_vcf.input.annotation,
        panel=rules.extract_panel_genes_from_vcf.input.panel,
    output:
        fasta=RESULTS / "drprg/popn_prg/genes.fa",
        faidx=RESULTS / "drprg/popn_prg/genes.fa.fai",
    log:
        LOGS / "create_references.log",
    params:
        padding=config.get("padding", 100),
    conda:
        str(ENVS / "create_references.yaml")
    script:
        str(SCRIPTS / "create_references.py")
