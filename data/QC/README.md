## Pipeline overview

![rulegraph](resources/rulegraph.png)

This directory houses the pipeline to perform quality control on the samples. There is
an assumption that the [Illumina
baseline variant calling](../../analysis/baseline_variants) has been done already - as
these VCFs are used for lineage-calling.

In the end of this pipeline the FASTQ files in the `subsampled/` directory are
decontaminated and contain no unmapped reads (as determined by the decontamination
database). These are the FASTQs used for all subsequent analyses.

## Results

The resulting report with sample composition and coverage is in
[`report.html.gz`](./report/report.html.gz). You will need to decompress the file before
viewing it. Something like `gzip -d -c report.html.gz > report.html` should do the
trick.

A description of why samples have been deemed to have failed QC can be found in
[`failed_qc_reasons.md`](../../docs/failed_qc_reasons.md).
