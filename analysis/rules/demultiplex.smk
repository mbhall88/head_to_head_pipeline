def determine_demultiplex_action(wildcards, input, output, threads, resources):
    df = samples.loc[(wildcards.region, wildcards.run)]
    out_dir = Path("analysis/{region}/nanopore/{run}/demultiplex/".format(region=wildcards.region, run=wildcards.run))

    expected_barcodes = df["nanopore_barcode"]
    is_multiplexed = any(expected_barcodes.isnull())

    result = {
        "df": df,
        "out_dir" = out_dir,
        "classification_path": out_dir / "classifications.txt",
        "is_multiplexed": is_multiplexed,
        "expected_barcodes": set(expected_barcodes)
    }
    return result


rule demultiplex:
    input:
        fast5s = "data/{region}/nanopore/{run}/f5s",
        fastq = "analysis/{region}/nanopore/{run}/basecalled.fastq"
    output:
        "analysis/{region}/nanopore/{run}/demultiplex/DEMULTIPLEX_COMPLETE"
    threads:
        config["demultiplex"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["demultiplex"]["memory"]
    singularity:
        config["demultiplex"]["container"]
    params:
        option = determine_demultiplex_action
    log:
        "analysis/logs/demultiplex_{region}_{run}.log"
    shell:
        """
        bash ../scripts/demultiplex.sh {params.option.is_multiplexed} \
            {input.fast5s} \
            {params.option.classification_path} \
            {input.fastq} \
            {params.option.out_dir} \
            {output} 2> {log}
        """


rule fix_filenames:
    input:
        "analysis/{region}/nanopore/{run}/demultiplex/DEMULTIPLEX_COMPLETE"
    output:
        "analysis/{region}/nanopore/{run}/demultiplex/FIX_NAMES_COMPLETE"
    threads: 1
    resources:
        mem_mb = 500
    params:
        option = determine_demultiplex_action
    log:
        "analysis/logs/fix_filenames_{region}_{run}.log"
    run:
        df = params.option.df
        out_dir = params.option.out_dir
        if params.option.is_multiplexed:
            # change names to sample IDs and remove redundant barcode
            expected_barcodes = params.option.expected_barcodes
            for fq in list(out_dir.rglob("*.fastq.gz")):
                barcode = fq.stem.split(".")[0]
                if barcode in expected_barcodes:
                    sample_id = df[df["nanopore_barcode"] == barcode].iloc[0]["sample_id"]
                    new_fname = out_dir / (sample_id + ".fastq.gz")
                    fq.replace(new_fname)
                else:
                    fq.unlink()  # delete file
        else:
            # copy basecalled fastq to demultiplexed folder with name of sample ID
            sample_id = df["sample_id"]
            if len(sample_id) > 1:
                raise ValueError(
                    "{} {} is not barcoded but has multiple sample IDs: {}".format(
                        wildcards.region, wildcards.run, sample_id
                    )
                )

            shell("gzip -c {input.fastq} > {out_dir}/{sample_id[0]}.fastq.gz 2> {log}")
        shell("touch {output}")
