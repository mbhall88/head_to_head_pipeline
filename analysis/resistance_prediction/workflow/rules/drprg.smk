rule download_panel:
    output:
        panel=RESULTS / "drprg/panel/panel.tsv",
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
