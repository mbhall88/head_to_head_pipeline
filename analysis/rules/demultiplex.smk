

rule demultiplex:
    input:
        fast5s = "data/{region}/nanopore/{run}/f5s",
        fastq = "analysis/{region}/nanopore/{run}/basecalled.fastq"
    output:
        "analysis/{region}/nanopore/{run}/demultiplex/{sample_id}.fastq.gz"
    params:
        out_dir = "analysis/{region}/nanopore/{run}/demultiplex/"
    log:
        "analysis/logs/demultiplex_{region}_{run}_{sample_id}.log"
    run:
        from pathlib import Path

        df = samples.loc[(wildcards.region, wildcards.run)]

        expected_barcodes = df["nanopore_barcode"]
        is_multiplexed = any(expected_barcodes.isnull())

        if is_multiplexed:
            classification_path = Path(params.out_dir) / "classifications.txt"
            # deepbinner classify
            shell(
                """deepbinner classify
                    --native {input.fast5s} > {classification_path} 2> {log}"""
                    )
            # deepbinner bin
            shell(
                """deepbinner bin --classes {classifications_path}
                    --reads {input.fastq}
                    --out_dir {params.out_dir} 2> {log}"""
                )
            # change names to sample IDs and remove redundant barcode
            expected_barcodes = set(expected_barcodes)
            for fq in list(Path(params.out_dir).rglob("*.fastq.gz")):
                barcode = fq.stem.split('.')[0]
                if barcode in expected_barcodes:
                    sample_id = df[df["nanopore_barcode"]==barcode].iloc[0]["sample_id"]
                    new_fname = Path(params.out_dir) / (sample_id + ".fastq.gz")
                    fq.replace(new_fname)
                else:
                    fq.unlink()  # delete file

        else:
            # copy basecalled fastq to demultiplexed folder with name of sample ID
            sample_id = df["sample_id"]
            if len(sample_id) > 1:
                raise ValueError("{} {} is not barcoded but has multiple sample IDs: {}".format(wildcards.region, wildcards.run, sample_id))

            shell("gzip -c {input.fastq} > {output} 2> {log}")



# add a rule that changes the names of the barcode files to their sample ID
# need to add some logic in there that will combine all the barcode files in
# non barcoded samples into a single file using the following logic
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rules
