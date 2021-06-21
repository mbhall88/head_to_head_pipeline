def infer_reads(wildcards, merged: bool = False):
    site = wildcards.site
    sample = wildcards.sample
    tech = wildcards.tech

    if tech == "nanopore":
        return f"subsampled/{site}/{tech}/{sample}/{sample}.subsampled.fastq.gz"
    elif merged: #todo
    else:
        return [
            f"subsampled/{site}/{tech}/{sample}/{sample}.subsampled.R{i}.fastq.gz"
            for i in [1, 2]
        ]

