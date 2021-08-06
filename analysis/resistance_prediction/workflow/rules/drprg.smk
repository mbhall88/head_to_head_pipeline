from pathlib import Path


rule download_panel:
    output:
        panel=RESULTS / "drprg/panel/panel.original.tsv",
    resources:
        mem_mb=int(0.5 * GB),
    log:
        LOGS / "download_panel.log",
    container:
        CONTAINERS["base"]
    params:
        url=config["mykrobe_panel_url"],
    shell:
        "wget -O {output.panel} {params.url} 2> {log}"


rule add_drugs_to_panel:
    input:
        panel=rules.download_panel.output.panel,
    output:
        panel=RESULTS / "drprg/panel/panel.resistant.tsv",
    resources:
        mem_mb=int(0.5 * GB),
    log:
        LOGS / "add_drugs_to_panel.log",
    container:
        CONTAINERS["conda"]
    params:
        mapping_url="https://ndownloader.figshare.com/files/14120195",
        delim=",",
    script:
        str(SCRIPTS / "add_drugs_to_panel.py")


rule add_known_non_resistance_vars_to_panel:
    input:
        panel=rules.add_drugs_to_panel.output.panel,
        known=RESOURCES / "non_resistant.tsv",
    output:
        panel=RESULTS / "drprg/panel/panel.tsv",
    log:
        LOGS / "add_known_non_resistance_vars_to_panel.log",
    container:
        CONTAINERS["base"]
    resources:
        mem_mb=int(0.3 * GB),
    shell:
        "awk 1 {input.panel} {input.known} > {output.panel} 2> {log}"


rule filter_panel:
    input:
        panel=rules.download_panel.output.panel,
    output:
        panel=RESULTS / "drprg/panel/panel.filtered.tsv",
    log:
        LOGS / "filter_panel.log",
    resources:
        mem_mb=int(0.5 * GB),
    params:
        exclude=["pncA", "katG"],
        frameshift_lengths=[1, 2],
    container:
        CONTAINERS["conda"]
    script:
        str(SCRIPTS / "filter_panel.py")


rule drprg_build:
    input:
        panel=rules.add_known_non_resistance_vars_to_panel.output.panel,
        ref=RESOURCES / "h37rv.fa",
        annotation=RESOURCES / "h37rv.gff3",
        prg_index=rules.index_popn_prg.output.index,
    output:
        outdir=directory(RESULTS / "drprg/index"),
        prg=RESULTS / "drprg/index/dr.prg",
        vcf=RESULTS / "drprg/index/panel.bcf",
        vcf_idx=RESULTS / "drprg/index/panel.bcf.csi",
        ref=RESULTS / "drprg/index/genes.fa",
    log:
        LOGS / "drprg_build.log",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    threads: 8
    container:
        CONTAINERS["drprg"]
    params:
        options=" ".join(
            [
                "-v",
                f"--match-len {config['match_len']}",
                f"--padding {config['padding']}",
            ],
        ),
        prebuilt_dir=lambda wildcards, input: Path(input.prg_index).parent,
    shell:
        """
        drprg build {params.options} -a {input.annotation} -o {output.outdir} -i {input.panel} \
          -f {input.ref} -t {threads}  -d {params.prebuilt_dir} 2> {log}
        """


rule seqtk_mergepe:
    input:
        reads=lambda wildcards: infer_reads(wildcards),
    output:
        merged=(
            RESULTS / "drprg/mergepe/{tech}/{site}/{sample}/{sample}.merged.fastq.gz"
        ),
    params:
        compress_lvl=9,
    log:
        LOGS / "seqtk_mergepe/{tech}/{site}/{sample}.log",
    threads: 2
    wrapper:
        "0.75.0-7-g05edf56/bio/seqtk/mergepe"


rule drprg_predict:
    input:
        index=rules.drprg_build.output.outdir,
        reads=lambda wildcards: infer_reads(wildcards, merged=True),
    output:
        outdir=directory(RESULTS / "drprg/predict/{tech}/{site}/{sample}"),
        report=RESULTS / "drprg/predict/{tech}/{site}/{sample}/{sample}.drprg.json",
        vcf=RESULTS / "drprg/predict/{tech}/{site}/{sample}/{sample}.drprg.bcf",
    log:
        LOGS / "drprg_predict/{tech}/{site}/{sample}.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    container:
        CONTAINERS["drprg"]
    params:
        opts=" ".join(["--verbose", "-s {sample}", "-u", "--failed"]),
        filters=lambda wildcards: drprg_filter_args(wildcards),
        tech_flag=lambda wildcards: "-I" if wildcards.tech == "illumina" else "",
    shell:
        """
        drprg predict {params.opts} {params.filters} {params.tech_flag} \
          -o {output.outdir} -i {input.reads} -x {input.index} -t {threads} 2> {log}
        """


rule subtract_panel_variants_from_drprg:
    input:
        panel=rules.drprg_build.output.vcf,
        query=rules.drprg_predict.output.vcf,
    output:
        vcf=RESULTS / "drprg/filtered_vcfs/{tech}/{site}/{sample}.novel.vcf",
    resources:
        mem_mb=int(0.5 * GB),
    log:
        LOGS / "subtract_panel_variants_from_drprg/{tech}/{site}/{sample}.log",
    conda:
        str(ENVS / "subtract_variants.yaml")
    script:
        str(SCRIPTS / "subtract_variants.py")
