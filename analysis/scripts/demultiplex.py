import subprocess
from pathlib import Path

params = snakemake.params
samples = params.samples
wildcards = snakemake.wildcards
log = snakemake.log

df = samples.loc[(wildcards.region, wildcards.run)]

expected_barcodes = df["nanopore_barcode"]
is_multiplexed = any(expected_barcodes.isnull())

if is_multiplexed:
    classification_path = Path(params.out_dir) / "classifications.txt"
    # deepbinner classify
    subprocess.run(
        "deepbinner classify --native {fast5s} > {classification_path} 2> {log}".format(
            fast5s=snakemake.input.fast5s,
            classification_path=classification_path,
            log=log,
        )
    )
    # deepbinner bin
    subprocess.run(
        "deepbinner bin --classes {classifications_path} --reads {fastq} --out_dir {out_dir} 2> {log}".format(
            classification_path=classification_path,
            fastq=snakemake.input.fastq,
            out_dir=params.out_dir,
            log=log,
        )
    )
    # change names to sample IDs and remove redundant barcode
    expected_barcodes = set(expected_barcodes)
    for fq in list(Path(params.out_dir).rglob("*.fastq.gz")):
        barcode = fq.stem.split(".")[0]
        if barcode in expected_barcodes:
            sample_id = df[df["nanopore_barcode"] == barcode].iloc[0]["sample_id"]
            new_fname = Path(params.out_dir) / (sample_id + ".fastq.gz")
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

    subprocess.run(
        "gzip -c {fastq} > {out_dir}/{sample_id}.fastq.gz 2> {log}".format(
            fastq=snakemake.input.fastq,
            outdir=params.out_dir,
            sample_id=sample_id[0],
            log=log,
        )
    )

subprocess.run("touch {output}".format(output=snakemake.output))
