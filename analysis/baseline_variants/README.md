The truth set of variants for the Illumina data will come from the PHE variant calling
pipeline `compass`.  
As a baseline for the Nanopore data, we will run `bcftools`, with some filtering of
variants to remove low-quality calls and with a mask to avoid repetitive and
structurally variable regions of the Mtb genome.

The results for this subsection will show concordance of the basic Nanopore variant
calling with the Illumina truth set. We will discuss how the different filters trialled
affected the results and plot recall against concordance.

## Pipeline overview

![rulegraph](resources/rulegraph.png)


## Results

The [report](./report.html) contains plots of call rate, concordance, and depth against
each other. In addition, it contains JSON files for each sample with it's concordance
and call rate raw values.
