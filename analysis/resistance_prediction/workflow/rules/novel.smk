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
        adjust_pos=True,
    conda:
        str(ENVS / "extract_panel_genes_from_vcf.yaml")
    script:
        str(SCRIPTS / "extract_panel_genes_from_vcf.py")


rule index_compass_panel_vcf:
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
        query_idx=rules.index_compass_panel_vcf.output[0],
        query=rules.extract_panel_genes_from_compass_vcf.output.vcf,
    output:
        vcf=RESULTS / "novel/filtered_vcfs/{sample}.novel.bcf",
    resources:
        mem_mb=int(0.5 * GB),
    log:
        LOGS / "subtract_panel_variants_from_compass/{sample}.log",
    conda:
        str(ENVS / "subtract_variants.yaml")
    script:
        str(SCRIPTS / "subtract_variants.py")


rule index_final_compass_panel_vcf:
    input:
        rules.subtract_panel_variants_from_compass.output.vcf,
    output:
        RESULTS / "novel/filtered_vcfs/{sample}.novel.bcf.csi",
    params:
        extra="-f",
    wrapper:
        "0.77.0/bio/bcftools/index"


rule index_drprg_ref_genes:
    input:
        rules.drprg_build.output.ref,
    output:
        RESULTS / "drprg/index/genes.fa.faidx",
    log:
        LOGS / "index_drprg_ref_genes.log",
    resources:
        mem_mb=int(0.3 * GB),
    params:
        "",
    wrapper:
        "0.77.0/bio/samtools/faidx"


rule assess_drprg_novel_calls:
    input:
        truth_vcf=rules.subtract_panel_variants_from_compass.output.vcf,
        truth_idx=rules.index_final_compass_panel_vcf.output[0],
        query_vcf=rules.subtract_panel_variants_from_drprg.output.vcf,
        query_idx=rules.index_novel_drprg_vcf.output[0],
        ref=rules.index_drprg_ref_genes.input[0],
        ref_idx=rules.index_drprg_ref_genes.output[0],
    output:
        summary=(
            RESULTS / "novel/assessment/{tech}/{site}/{sample}/{sample}.summary.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: int(GB) * attempt,
    log:
        LOGS / "assess_drprg_novel_calls/{tech}/{site}/{sample}.log",
    container:
        CONTAINERS["happy"]
    params:
        opts=" ".join(
            (
                "--set-gt hom",
                "--pass-only",
                "--write-vcf",
                "--leftshift",
            )
        ),
        prefix=lambda wc, output: output.summary.split(".")[0],
    shell:
        """
        truth_count=$(bcftools view -f -H {input.truth_vcf} | wc -l)
        query_count=$(bcftools view -f -H {input.query_vcf} | wc -l)

        if [ "$truth_count" -eq 0 ] || [ "$query_count" -eq 0 ]; then
          printf 'TP,FN,FP\n0,%d,%d\n' "$truth_count" "$query_count" > {output.summary} 2> {log}
        else
          hap.py {params.opts} -o {params.prefix} -r {input.ref} \
            {input.truth_vcf} {input.query_vcf} 2> {log}
        fi
        """
