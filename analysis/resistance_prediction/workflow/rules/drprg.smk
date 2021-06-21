rule download_panel:
    output:
        panel=RESULTS / "drprg/panel/panel.original.tsv",
    resources:
        mem_mb=int(0.5 * GB),
    log:
        LOGS / "download_panel.log",
    params:
        url=config["mykrobe_panel_url"],
    shell:
        "wget -O {output.panel} {params.url} 2> {log}"


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


rule add_drugs_to_panel:
    input:
        panel=rules.filter_panel.output.panel,
    output:
        panel=RESULTS / "drprg/panel/panel.tsv",
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


rule drprg_build:
    input:
        panel=rules.add_drugs_to_panel.output.panel,
        ref=RESOURCES / "h37rv.fa",
        annotation=RESOURCES / "h37rv.gff3",
    output:
        outdir=directory(RESULTS / "drprg/index"),
        prg=RESULTS / "drprg/index/dr.prg",
    log:
        LOGS / "drprg_build.log",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    threads: 8
    container:
        CONTAINERS["drprg"]
    params:
        " ".join(
            [
                "-v",
                "--force",
                f"--match-len {config['match_len']}",
                f"--padding {config['padding']}",
            ],
        ),
    shell:
        """
        drprg build {params} -a {input.annotation} -o {output.outdir} -i {input.panel} \
          -f {input.ref} -t {threads} 2> {log}
        """


rule drprg_predict:
    input:
        index=rules.drprg_build.output.outdir,
        reads=lambda wildcards: QC(infer_reads(wildcards)),
    output:
        outdir=directory(RESULTS / "drprg/predict/{tech}/{site}/{sample}"),
        report=RESULTS / "drprg/predict/{tech}/{site}/{sample}/{sample}.drprg.json",
    log:
        LOGS / "drprg_predict/{tech}/{site}/{sample}.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    container:
        CONTAINERS["drprg"]
    params:
        opts=" ".join(
            [
                "--verbose",
                "--force",
                "-s {sample}",
                f"-d {config['filters']['min_covg']}",
                f"-b {config['filters']['min_strand_bias']}",
                f"-g {config['filters']['min_gt_conf']}",
                f"-L {config['filters']['max_indel']}",
                f"-K {config['filters']['min_frs']}",
            ]
        ),
        tech_flag=lambda wildcards: "-I" if wildcards.tech == "illumina" else "",
    shell:
        """
        drprg predict {params.opts} {params.tech_flag} -o {output.outdir} \
          -i {input.reads} -x {input.index} -t {threads} 2> {log}
        """
