def infer_reads(wildcards):
    site = wildcards.site
    sample = wildcards.sample
    tech = wildcards.tech

    if tech == "nanopore":
        return f"subsampled/{site}/{tech}/{sample}/{sample}.subsampled.fastq.gz"
    else:
        return [
            f"subsampled/{site}/{tech}/{sample}/{sample}.subsampled.R{i}.fastq.gz"
            for i in [1, 2]
        ]

