def infer_reads(wildcards, merged: bool = False):
    site = wildcards.site
    sample = wildcards.sample
    tech = wildcards.tech

    if tech == "nanopore":
        return QC(f"subsampled/{site}/{tech}/{sample}/{sample}.subsampled.fastq.gz")
    elif merged:
        return (
            RESULTS / f"drprg/mergepe/{tech}/{site}/{sample}/{sample}.merged.fastq.gz"
        )
    else:
        return QC([
            f"subsampled/{site}/{tech}/{sample}/{sample}.subsampled.R{i}.fastq.gz"
            for i in [1, 2]
        ])
