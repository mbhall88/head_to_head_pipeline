# *Mycobacterium tuberculosis* population reference graph construction

todo: add rulegraph

## Introduction

The code in this section of the repository is for construction of the
"basic" H37Rv PRG.

We refer to this as "basic" because we will be using H37Rv as the
skeleton, and it is known that it misses some genes found in other
lineages.

The outline for the pipeline will be:
- Split H37Rv into chunks based on the annotation.
- Construct merged VCFs for each lineage for all VCFs in the CRyPTIC
  dataset (17/03/2020)
- For a range of allele frequencies (AFs), apply variants from each VCF
  to their respective H37Rv chunks. By apply, we mean, make a copy of
  the reference sequence for that chunk, apply the variant to that copy
  and add it to the chunk to make a multi-sequence FASTA file.
- Run [`make_prg`][make_prg] on all chunks.
- Combine resulting PRGs into a single PRG file for use with
  [`pandora`][pandora].


[make_prg]: https://github.com/rmcolq/make_prg
[pandora]: https://github.com/rmcolq/pandora

