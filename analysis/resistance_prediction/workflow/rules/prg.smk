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


rule index_popn_vcf:
    input:
        rules.extract_panel_genes_from_vcf.output.vcf,
    output:
        RESULTS / "drprg/popn_prg/popn.bcf.csi",
    params:
        extra="-f",
    wrapper:
        "0.76.0/bio/bcftools/index"


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


rule create_popn_pre_msas:
    input:
        vcf=rules.extract_panel_genes_from_vcf.output,
        vcfidx=rules.index_popn_vcf.output[0],
        references=rules.create_references.output.fasta,
    output:
        directory(RESULTS / "drprg/popn_prg/pre_msas"),
    log:
        LOGS / "create_popn_pre_msas.log",
    conda:
        str(ENVS / "create_popn_pre_msas.yaml")
    script:
        str(SCRIPTS / "create_popn_pre_msas.py")


rule create_popn_msas:
    input:
        pre_msas=rules.create_popn_pre_msas.output[0],
    output:
        directory(RESULTS / "drprg/popn_prg/msas"),
    log:
        LOGS / "create_popn_msas.log",
    container:
        CONTAINERS["mafft"]
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
    shell:
        """
        mkdir -p {output[0]} 2> {log}
        for f in {input.pre_msas}/*.fa
        do
            outname={output[0]}/$(basename "$f")
            mafft --auto --thread {threads} "$f" > "$outname"
        done 2>> {log}
        """
