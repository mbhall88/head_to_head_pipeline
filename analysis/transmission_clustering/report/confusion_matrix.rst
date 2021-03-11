The confusion matrix from comparing whether two samples (within an Illumina SNP distance of {{ snakemake.wildcards.threshold }})
are clustered together by both COMPASS (Illumina) and {{ snakemake.wildcards.caller }} (Nanopore). We classify each pair
according to the following definition:

- **True Positive (TP)**: the pair **is** clustered in the truth (COMPASS) and also in the predicted (i.e. {{ snakemake.wildcards.caller }})
- **True Negative (TN)**:  the pair **is not** clustered in the truth and also in the predicted
- **False Positive (FP)**: the pair **is not** clustered in the truth but **is** clustered in the predicted
- **False Negative (FN)**: the pair **is** clustered in the truth but **is not** clustered in the predicted