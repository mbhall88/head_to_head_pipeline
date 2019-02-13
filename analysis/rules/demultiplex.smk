rule demultiplex:
    input:
        fast5s = "data/{region}/nanopore/{run}/f5s",
        fastq = "analysis/{region}/nanopore/{run}/basecalled.fastq"
    output:
        "analysis/{region}/nanopore/{run}/demultiplex/COMPLETE"
    threads:
        config["demultiplex"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["demultiplex"]["memory"]
    singularity:
        config["demultiplex"]["container"]
    params:
        out_dir = "analysis/{region}/nanopore/{run}/demultiplex/",
        samples = samples
    log:
        "analysis/logs/demultiplex_{region}_{run}.log"
    script:
        "analysis/scripts/demultiplex.py"





# add a rule that changes the names of the barcode files to their sample ID
# need to add some logic in there that will combine all the barcode files in
# non barcoded samples into a single file using the following logic
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rules
