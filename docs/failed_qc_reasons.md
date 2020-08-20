The column 'failed_qc' in [`samplesheet.csv`](samplesheet.csv) indicates samples that
have a `1` in the 'excluded' column, or samples that failed [QC](data/QC). The reasons
for exclusion prior to QC can be found in
[`exclusion_reasons.md`](exclusion_reasons.md). This document is for describing why
samples have failed QC.


[TOC]: #

# Table of Contents
- [R23399](#r23399)
- [R28773](#r28773)
- [R21046](#r21046)
- [mada_1-9](#mada_1-9)
- [mada_101](#mada_101)
- [mada_153](#mada_153)
- [mada_145](#mada_145)
- [mada_1-52](#mada_1-52)
- [mada_1-26](#mada_1-26)
- [mada_1-23](#mada_1-23)
- [mada_1-29](#mada_1-29)
- [mada_1-49](#mada_1-49)
- [R20968](#r20968)
- [R22333](#r22333)
- [R22727](#r22727)
- [R21648](#r21648)
- [R37259](#r37259)
- [R20328](#r20328)
- [R16657](#r16657)
- [R36440](#r36440)
- [R24409](#r24409)
- [R27273](#r27273)
- [R29669](#r29669)
- [R22393](#r22393)
- [R24644](#r24644)
- [R17830](#r17830)
- [R28195](#r28195)
- [R26792](#r26792)
- [R23144](#r23144)
- [R21911](#r21911)
- [R28538](#r28538)
- [R20412](#r20412)
- [R28691](#r28691)
- [R38290](#r38290)
- [mada_114](#mada_114)
- [mada_146](#mada_146)
- [mada_1-21](#mada_1-21)
- [mada_2-47](#mada_2-47)
- [mada_1-42](#mada_1-42)
- [mada_1-37](#mada_1-37)
- [mada_1-31](#mada_1-31)
- [mada_1-35](#mada_1-35)
- [mada_1-24](#mada_1-24)
- [mada_147](#mada_147)
- [18_0622317](#18_0622317)
- [18_0622293](#18_0622293)
- [18_0622212](#18_0622212)
- [18_0622034](#18_0622034)
- [18_620602](#18_620602)
- [17_616166](#17_616166)
- [18_0622441](#18_0622441)
- [18_620650](#18_620650)
- [18_0622429](#18_0622429)
- [18_0621850](#18_0621850)
- [18_0622302](#18_0622302)
- [18_0621852](#18_0621852)
- [18_0622217](#18_0622217)


## R23399

No lineage could be assigned for this sample.

At the end of the QC pipeline, this sample only has 0.3x Illumina and 0.4x Nanopore
data. See https://github.com/mbhall88/head_to_head_pipeline/issues/23

*Note: This sample also had a very high fraction of reads (>85%) for both technologies
classified as unmapped when mapping to the decontamination database.*

## R28773

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

Nanopore coverage 10.15 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R21046

No lineage could be assigned for this sample.

*Note: this sample has good coverage for both technologies and very little
contamination. The logs for the lineage assignment indicate that at all
lineage-definining positions, the sample failed the FILTER with `A25;B1;K0.90;f0.35;n5`,
which effectively means it failed multiple quality requirements at each of these
positions.*

## mada_1-9

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 4.1.2
- 2.2.5
- Bovis
- 4.2
- 2.2.7
- 4.4.1.2
- 4.6.2.1
- 3.1.2
- 4.1.1.3

Additionally, after filtering and subsampling, the Illumina data only had 12.1x
coverage.

*Note: more than half the reads for both technologies mapped to contamination and more
reads were classified as unmapped than to keep.*

## mada_101

No lineage could be assigned for this sample.

*Note: this samples has a lot of coverage for both technologies and very little
contamination. The logs for the lineage assignment indicate it had a few of the
lineage-defining positions fail the COMPASS filter where less than 90% of basecalls were
in agreement at this position. As such those positions were het calls.*

## mada_153

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 4
- 3

*Note: This sample has extremely high coverage for both technologies with very little
contamination. Looking at the lineage assignment logs, it seems to have failed FILTER at
two positions due to low support.*

## mada_145

No lineage could be assigned for this sample.

*Note: this samples has a lot of coverage for both technologies and very little
contamination. The logs for the lineage assignment indicate it had a number of the
lineage-defining positions fail the COMPASS filter where less than 90% of basecalls were
in agreement at this position. As such those positions were het calls.*

## mada_1-52

No lineage could be assigned for this sample.

*Note: this samples has a lot of coverage for both technologies and very little
contamination. The logs for the lineage assignment indicate it had a number of the
lineage-defining positions fail the COMPASS filter where less than 90% of basecalls were
in agreement at this position. As such those positions were het calls.*

## mada_1-26

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 2.2.5
- 4.2
- 2.2.7
- 4.6.2.1
- 4.6.1

At the end of the QC pipeline, this sample only has 2x Illumina and 2.8x Nanopore data.
See https://github.com/mbhall88/head_to_head_pipeline/issues/23

*Note: more than half the reads for both technologies mapped to contamination and more
reads were classified as unmapped than to keep.*

## mada_1-23

At the end of the QC pipeline, this sample only has 0.2x Illumina and 0.6x Nanopore
data. See https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_1-29

At the end of the QC pipeline, this sample only has 0.06x Illumina and 0.46x Nanopore
data. See https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_1-49

At the end of the QC pipeline, this sample only has 0.07x Illumina and 1.83x Nanopore
data. See https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R20968

Nanopore coverage 24.87 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R22333

Nanopore coverage 1.21 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R22727

Nanopore coverage 12.19 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R21648

Nanopore coverage 5.28 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R37259

Nanopore coverage 22.44 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R20328

Nanopore coverage 0.41 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R16657

Nanopore coverage 27.62 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R36440

Nanopore coverage 18.33 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R24409

Nanopore coverage 9.75 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R27273

Nanopore coverage 22.86 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R29669

Nanopore coverage 26.02 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R22393

Nanopore coverage 20.73 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R24644

Nanopore coverage 19.86 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R17830

Nanopore coverage 0.01 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R28195

Nanopore coverage 22.32 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R26792

Nanopore coverage 4.56 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R23144

Nanopore coverage 20.53 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R21911

Nanopore coverage 6.3 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R28538

Nanopore coverage 0.0 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R20412

Nanopore coverage 4.32 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R28691

Nanopore coverage 21.97 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## R38290

Nanopore coverage 4.17 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_114

Nanopore coverage 1.01 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_146

Nanopore coverage 0.67 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_1-21

Nanopore coverage 29.77 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_2-47

Nanopore coverage 21.9 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_1-42

Nanopore coverage 6.91 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_1-37

Nanopore coverage 10.74 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_1-31

Nanopore coverage 15.7 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_1-35

Nanopore coverage 9.42 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_1-24

Nanopore coverage 15.94 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## mada_147

Nanopore coverage 0.28 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## 18_0622317

Nanopore coverage 22.29 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## 18_0622293

Nanopore coverage 22.04 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## 18_0622212

Nanopore coverage 8.97 is below cutoff of 30x. See
<https://github.com/mbhall88/head_to_head_pipeline/issues/23>. In addition, the Illumina
data for this sample, 16.28x, is below the threshold of 20x stated in
<https://github.com/mbhall88/head_to_head_pipeline/issues/23#issuecomment-660900615>.

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

## 18_0622034

Nanopore coverage 25.13 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## 18_620602

Nanopore coverage 22.42 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

## 17_616166

Nanopore coverage 11.08 is below cutoff of 30x. See
https://github.com/mbhall88/head_to_head_pipeline/issues/23

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

## 18_0622441

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

## 18_620650

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

*Note: This sample has extremely high coverage for both technologies with very little
contamination.*

## 18_0622429

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

*Note: This sample has extremely high coverage for both technologies with very little
contamination.*

## 18_0621850

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

*Note: This sample has extremely high coverage for both technologies with very little
contamination.*

## 18_0622302

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

*Note: This sample has extremely high coverage for both technologies with very little
contamination.*

## 18_0621852

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

*Note: This sample has extremely high coverage for both technologies with very little
contamination.*

## 18_0622217

This sample's lineage assignment was considered `mixed`. The found lineages were:
- 3
- 4

*Note: This sample has extremely high coverage for both technologies but shows about 5%
contamination in both technologies.*
