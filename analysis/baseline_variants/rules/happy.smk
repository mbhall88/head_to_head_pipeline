rule compress_and_index_truth_vcf:
    input:
        vcf=rules.evaluate_compass_snps_with_truth_assembly.output.truth_vcf,
    output:
        bcf=truth_eval_dir / "varifier/{sample}/compass/recall/truth_vcf/04.truth.bcf",
        index=truth_eval_dir / "varifier/{sample}/compass/recall/truth_vcf/04.truth.bcf.csi",
    container:
        containers["bcftools"]
    log:
        rule_log_dir / "compress_truth_vcf/{sample}.log"
    shell:
        """
        bcftools view -O b -o {output.bcf} {input.vcf} 2> {log}
        bcftools index -f -o {output.index} {output.bcf} 2>> {log}
        """


rule index_compass_vcf:
    input:
        vcf=compass_vcf_dir / "{sample}.compass.vcf.gz",
    output:
        index=compass_vcf_dir / "{sample}.compass.vcf.gz.csi",
    container:
        containers["bcftools"]
    log:
        rule_log_dir / "index_compass_vcf/{sample}.log"
    shell:
        "bcftools index -f -o {output.index} {input.vcf} 2> {log}"


rule index_bcftools_vcf:
    input:
        vcf=nanopore_dir / "filtered_snps/madagascar/{sample}.snps.filtered.bcf",
    output:
        index=nanopore_dir / "filtered_snps/madagascar/{sample}.snps.filtered.bcf.csi",
    container:
        containers["bcftools"]
    log:
        rule_log_dir / "index_bcftools_vcf/{sample}.log"
    shell:
        "bcftools index -f -o {output.index} {input.vcf} 2> {log}"


rule index_truth_assembly:
    input:
        asm=asm_dir / "{sample}/flye/pacbio/decontam.assembly.flye.pacbio.fasta",
    output:
        index=asm_dir / "{sample}/flye/pacbio/decontam.assembly.flye.pacbio.fasta.fai",
    container:
        containers["samtools"]
    log:
        rule_log_dir / "index_truth_assembly/{sample}.log"
    shell:
        "samtools faidx {input.asm} 2> {log}"


rule evaluate_compass_with_happy:
    input:
        truth_vcf=rules.compress_and_index_truth_vcf.output.bcf,
        truth_idx=rules.compress_and_index_truth_vcf.output.index,
        query_vcf=rules.index_compass_vcf.input.vcf,
        query_idx=rules.index_compass_vcf.output.index,
        ref=rules.index_truth_assembly.input.asm,
        ref_idx=rules.index_truth_assembly.output.index,
        mask=asm_dir / "{sample}/flye/pacbio/assessment/{sample}.flye.accuracy.pacbio.bed",
    output:
        summary=(
            truth_eval_dir / "happy/{sample}/compass/{sample}.summary.csv",
        ),
    resources:
        mem_mb=lambda wildcards, attempt: int(8 * GB) * attempt,
    threads:
        8
    log:
        rule_log_dir / "evaluate_compass_with_happy/{sample}.log",
    container:
        containers["happy"]
    params:
        opts=" ".join(
            (
                # "--set-gt hom",
                "--pass-only",
                "--write-vcf",
                "--leftshift",
                "--engine=vcfeval"
            )
        ),
        prefix=lambda wildcards, output: str(output.summary).split(".")[0],
    shell:
        """
        hap.py {params.opts} -o {params.prefix} --threads {threads} -r {input.ref} \
          -T ^{input.mask} {input.truth_vcf} {input.query_vcf} 2> {log}
        """

rule evaluate_bcftools_with_happy:
    input:
        truth_vcf=rules.compress_and_index_truth_vcf.output.bcf,
        truth_idx=rules.compress_and_index_truth_vcf.output.index,
        query_vcf=rules.index_bcftools_vcf.input.vcf,
        query_idx=rules.index_bcftools_vcf.output.index,
        ref=rules.index_truth_assembly.input.asm,
        ref_idx=rules.index_truth_assembly.output.index,
        mask=asm_dir / "{sample}/flye/pacbio/assessment/{sample}.flye.accuracy.pacbio.bed",
    output:
        summary=(
            truth_eval_dir / "happy/{sample}/bcftools/{sample}.summary.csv",
        ),
    resources:
        mem_mb=lambda wildcards, attempt: int(8 * GB) * attempt,
    threads:
        8
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
                "--engine=vcfeval"
            )
        ),
        prefix=lambda wildcards, output: str(output.summary).split(".")[0],
    shell:
        """
        hap.py {params.opts} -o {params.prefix} --threads {threads} -r {input.ref} \
          -T ^{input.mask} {input.truth_vcf} {input.query_vcf} 2> {log}
        """
