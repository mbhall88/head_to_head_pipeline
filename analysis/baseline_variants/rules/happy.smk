rule compress_and_index_truth_vcf:
    input:
        vcf=rules.evaluate_compass_snps_with_truth_assembly.output.truth_vcf,
    output:
        bcf=truth_eval_dir / "varifier/{sample}/compass/recall/truth_vcf/04.truth.bcf",
        index=truth_eval_dir
        / "varifier/{sample}/compass/recall/truth_vcf/04.truth.bcf.csi",
    container:
        containers["bcftools"]
    log:
        rule_log_dir / "compress_truth_vcf/{sample}.log",
    shell:
        """
        bcftools view -O b -o {output.bcf} {input.vcf} 2> {log}
        bcftools index -f -o {output.index} {output.bcf} 2>> {log}
        """


rule convert_and_index_compass_vcf:
    input:
        vcf=compass_vcf_dir / "{sample}.compass.vcf.gz",
    output:
        bcf=compass_vcf_dir / "{sample}.compass.bcf",
        index=compass_vcf_dir / "{sample}.compass.bcf.csi",
    container:
        containers["bcftools"]
    log:
        rule_log_dir / "index_compass_vcf/{sample}.log",
    shell:
        """
        bcftools view -O b -o {output.bcf} {input.vcf} 2> {log}
        bcftools index -f -o {output.index} {output.bcf} 2>> {log}
        """


rule index_bcftools_vcf:
    input:
        vcf=nanopore_dir / "filtered_snps/madagascar/{sample}.snps.filtered.bcf",
    output:
        index=nanopore_dir / "filtered_snps/madagascar/{sample}.snps.filtered.bcf.csi",
    container:
        containers["bcftools"]
    log:
        rule_log_dir / "index_bcftools_vcf/{sample}.log",
    shell:
        "bcftools index -f -o {output.index} {input.vcf} 2> {log}"


rule index_ref_genome:
    input:
        asm=H37RV["genome"],
    output:
        index=f"{H37RV['genome']}.fai",
    container:
        containers["samtools"]
    log:
        rule_log_dir / "index_truth_assembly.log",
    shell:
        "samtools faidx {input.asm} 2> {log}"


rule evaluate_compass_with_happy:
    input:
        truth_vcf=rules.compress_and_index_truth_vcf.output.bcf,
        truth_idx=rules.compress_and_index_truth_vcf.output.index,
        query_vcf=rules.convert_and_index_compass_vcf.output.bcf,
        query_idx=rules.convert_and_index_compass_vcf.output.index,
        ref=rules.index_ref_genome.input.asm,
        ref_idx=rules.index_ref_genome.output.index,
        mask=H37RV["mask"],
    output:
        summary=(truth_eval_dir / "happy/{sample}/compass/{sample}.summary.csv",),
        outdir=directory(truth_eval_dir / "happy/{sample}/compass"),
    resources:
        mem_mb=lambda wildcards, attempt: int(8 * GB) * attempt,
    threads: 8
    log:
        rule_log_dir / "evaluate_compass_with_happy/{sample}.log",
    container:
        containers["happy"]
    params:
        opts=" ".join(("--pass-only", "--write-vcf", "--leftshift")),
    shell:
        """
        hap.py {params.opts} -o {output.outdir}/{wildcards.sample} --threads {threads} \
          -r {input.ref} -T ^{input.mask} {input.truth_vcf} {input.query_vcf} 2> {log}
        """


rule evaluate_bcftools_with_happy:
    input:
        truth_vcf=rules.compress_and_index_truth_vcf.output.bcf,
        truth_idx=rules.compress_and_index_truth_vcf.output.index,
        query_vcf=rules.index_bcftools_vcf.input.vcf,
        query_idx=rules.index_bcftools_vcf.output.index,
        ref=rules.index_ref_genome.input.asm,
        ref_idx=rules.index_ref_genome.output.index,
        mask=H37RV["mask"],
    output:
        summary=(truth_eval_dir / "happy/{sample}/bcftools/{sample}.summary.csv",),
        outdir=directory(truth_eval_dir / "happy/{sample}/bcftools"),
    resources:
        mem_mb=lambda wildcards, attempt: int(8 * GB) * attempt,
    threads: 8
    log:
        rule_log_dir / "evaluate_bcftools_with_happy/{sample}.log",
    container:
        containers["happy"]
    params:
        opts=" ".join(
            (
                "--set-gt hom",
                "--pass-only",
                "--write-vcf",
                "--leftshift",
            )
        ),
    shell:
        """
        hap.py {params.opts} -o {output.outdir}/{wildcards.sample} --threads {threads} \
          -r {input.ref} -T ^{input.mask} {input.truth_vcf} {input.query_vcf} 2> {log}
        """
