The truth set of variants for the Illumina data will come from the PHE variant calling
pipeline `compass`.  
As a baseline for the Nanopore data, we will run `bcftools`, with some filtering of
variants to remove low-quality calls and with a mask to avoid repetitive and
structurally variable regions of the Mtb genome.

[TOC]: #

# Table of Contents
- [Pipeline overview](#pipeline-overview)
- [Methods](#methods)
  - [Illumina variant calls](#illumina-variant-calls)
  - [Nanopore variant calls](#nanopore-variant-calls)
  - [Filter Nanopore SNPs](#filter-nanopore-snps)
  - [Concordance of Nanopore calls to Illumina](#concordance-of-nanopore-calls-to-illumina)
  - [Truth evaluation with `varifier`](#truth-evaluation-with-varifier)
  - [Comparing SNP distances from Illumina and Nanopore](#comparing-snp-distances-from-illumina-and-nanopore)
- [Results](#results)
  - [Concordance](#concordance)
  - [Truth evaluation](#truth-evaluation)
  - [Distance](#distance)


## Pipeline overview

![rulegraph](resources/rulegraph.png)

## Methods

### Illumina variant calls

Variant calls for the Illumina data were done using the Public Health England pipeline
[compass]. This step is not included in the pipeline as it was run by Fan from Oxford.
See #17 and #36 for more information.

### Nanopore variant calls

Nanopore variants (SNPs) were called using [`bcftools`][bcftools]. Specifically, the
`mpileup` subcommand and then passing the output from that into the `call` subcommand.
Prior to `mpileup` the Nanopore reads were mapped to H37Rv with minimap2. The SAM files
produced from this mapping were the input to `mpileup`. When using the `call`
subcommand, we used the multiallelic caller<sup>[1]</sup>. Paramters used for this step
are linked:
[`mpileup`](https://github.com/mbhall88/head_to_head_pipeline/blob/cb61f1bd0b9db15d123c562f3135222ca9f9aba9/analysis/baseline_variants/Snakefile#L127-L148)
and
[`call`](https://github.com/mbhall88/head_to_head_pipeline/blob/cb61f1bd0b9db15d123c562f3135222ca9f9aba9/analysis/baseline_variants/Snakefile#L151-L167).

### Filter Nanopore SNPs

`bcftools` SNP calls were filtered using
[`apply_filters.py`](./scripts/apply_filters.py). The parameters for this step can be
found
[here](https://github.com/mbhall88/head_to_head_pipeline/blob/cb61f1bd0b9db15d123c562f3135222ca9f9aba9/analysis/baseline_variants/Snakefile#L170-L207).
The choice of filters was tightly linked to evaluation of the results from the
[truth evaluation](#truth-evaluation). A detailed view of how the final filters were
chosen can be seen in #39. The output from this step is a BCF file with the FILTER
column filled in with either PASS or the filter(s) that variant failed.

### Concordance of Nanopore calls to Illumina

How closely do the Nanopore variant calls mirror the Illumina calls from compass?

Definition of concordance-related terms (taken from
https://github.com/mbhall88/head_to_head_pipeline/issues/38#issuecomment-661780024) -
*note "true" is considered compass for this section of the analysis*:

**Call rate (for true ALTs)**

what % of true ALTs does bcftools make a REF/ALT call - anything except NULL and FILTER

**Concordance (for true ALTs)**

what % of true ALTs does bcftools genotype agree with truth (excludes NULL, FILTER, and
HET)

**Genomewide Call rate**

what % of true REF or ALT positions does bcftools make a REF/ALT call - anything except
NULL and FILTER

**Genomewide Concordance**

what % of true REF or ALT positions does bcftools genotype agree with truth (excludes
NULL, FILTER, and HET)

The configuration used to run this step of the pipeline can be found
[here](https://github.com/mbhall88/head_to_head_pipeline/blob/cb61f1bd0b9db15d123c562f3135222ca9f9aba9/analysis/baseline_variants/Snakefile#L210-L243).
The script used for calculating concordance is
[`concordance.py`](./scripts/concordance.py). The output from this step is, for each
sample, a CSV file containing a classification for each position contained in the union
of positions from both the Illumina and Nanopore VCFs. There is a classification for
each technology (i.e. REF, ALT, missing etc.) and a concordance classification, such as
TRUE_ALT, masked, FALSE_REF etc.

### Truth evaluation with `varifier`

An unbiased view of how correct the variant calls for each sequencing technology.

For the samples with good PacBio data, [`varifier`][varifier] was used to evaluate the
variant calls for both the Illumina and Nanopore data. These PacBio CCS assemblies are
considered "truth" and so will provide the best feedback method for fine-tuning the
Nanopore filters, and also for highlight the differences/similarities in variant calls
from Illumina and Nanopore reads. The parameters used for this step can be found
[here](https://github.com/mbhall88/head_to_head_pipeline/blob/cb61f1bd0b9db15d123c562f3135222ca9f9aba9/analysis/baseline_variants/Snakefile#L317-L430)
with further information about this step in #42 and #39.

### Comparing SNP distances from Illumina and Nanopore

Do Nanopore (`bcftools`) SNP calls provide similar pairwise distance values to Illumina
SNP calls?

The first step here was to generate a consensus sequence for each VCF file using
[`consensus.py`](./scripts/consensus.py).
[This rule](https://github.com/mbhall88/head_to_head_pipeline/blob/cb61f1bd0b9db15d123c562f3135222ca9f9aba9/analysis/baseline_variants/Snakefile#L443-L473)
outlines the parameters etc. used to run this section. Discussion of the choice of
parameters can be found in #41.  
When generating the consensus sequence, we effectively "mask" (replace reference base
with 'N') any positions that is in the mask, fails filters, is heterozygous, missing in
VCF, or is a null call.  
The pairwise SNP distance matrix was produced using [`snp-dists`][snp-dists] - see
[here](https://github.com/mbhall88/head_to_head_pipeline/blob/cb61f1bd0b9db15d123c562f3135222ca9f9aba9/analysis/baseline_variants/Snakefile#L497-L519)
for the usage.

## Results

The [report] contains results for the following parts of this analysis.

### Concordance

In this section of the [report], there is plots showing the relationship between call
rate and concordance and also depth with both call rate and concordance. For definitions
of these terms, see the [methods](#concordance-of-nanopore-calls-to-illumina).
Additionally, the JSON file summarising the concordance is attached to the report for
each sample.

### Truth evaluation

In this section of the [report], there is a plot of the precision and recall of the
compass and `bcftools` SNP calls. The raw summary statistics JSON files are also
included for each sample.

### Distance

This section of the [report] contains a heatmap for both compass and `bcftools` distance
matrices. It also contains two dotplots. Each dotplot shows the pairwise SNP distance
from compass (x-axis) against the distance for the same pair from `bcftools` (y-axis). A
linear equation is provided in the plots to describe the relationship. One dotplot is
for the full set of pairs, whilst the "close" dotplot is for pairs that have a (compass)
SNP distance no greater than 100.

[report]: ./report.html
[compass]: https://github.com/oxfordmmm/CompassCompact
[bcftools]: https://samtools.github.io/bcftools/bcftools.html
[varifier]: https://github.com/iqbal-lab-org/varifier
[snp-dists]: https://github.com/tseemann/snp-dists
[1]: https://samtools.github.io/bcftools/call-m.pdf

