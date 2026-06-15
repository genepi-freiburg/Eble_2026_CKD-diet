# Methods — MRAgent Validation Analysis

## Overview

We used MRAgent v0.2.5 (Xu et al., 2024; https://github.com/xuwei1997/MRAgent)
in its OE (Exposure–Outcome) mode to independently replicate all 10 × 2
exposure–outcome pairs from the primary analysis. The replication used an
MRAgent-compatible R script (`run_mragent_validation.R`) that implements
MRAgent's exact statistical workflow via `ieugwasr` directly.

## GWAS Sources

### Exposures (identical to primary analysis)

All exposure GWAS are UK Biobank dietary phenotypes, accessed via the IEU
OpenGWAS API (https://gwas.mrcieu.ac.uk/):

| Exposure | OpenGWAS ID | n |
|---|---|---|
| Tea intake | ukb-b-6066 | ~461 k |
| Salad/raw vegetable intake | ukb-b-1996 | ~461 k |
| Processed meat intake | ukb-b-6324 | ~461 k |
| Fresh fruit intake | ukb-b-3881 | ~461 k |
| Cooked vegetable intake | ukb-b-8089 | ~461 k |
| Dried fruit intake | ukb-b-16576 | ~461 k |
| Oily fish intake | ukb-b-2209 | ~461 k |
| Beef intake | ukb-b-2862 | ~461 k |
| Pork intake | ukb-b-5640 | ~461 k |
| Lamb/mutton intake | ukb-b-14179 | ~461 k |

### Outcomes

The primary analysis used local summary statistics downloaded directly from
CKDGen (https://ckdgen.imbi.uni-freiburg.de/). These are not fully available on
OpenGWAS (see note below). The closest available equivalents were used:

| Outcome | OpenGWAS ID | GWAS | n | Phenotype |
|---|---|---|---|---|
| eGFR | ebi-a-GCST90103634 | CKDGen 2019, European | 1,004,040 | eGFR (creatinine-based, log-transformed) |
| CKD | ebi-a-GCST90018822 | Sakaue et al. 2021 | 482,858 | Chronic renal failure (ICD-coded) |

**Note on CKD GWAS availability:**
The primary analysis used the CKDGen 2019 binary CKD GWAS (`ckd_EA.txt.gz`,
Wuttke et al. 2019, PMID 31152163). On OpenGWAS:
- `ieu-b-7440` (the nominal CKDGen 2019 CKD ID) returns empty results (API
  error — confirmed during analysis on 2026-05-05).
- `ieu-a-1102` (Pattaro 2015) covers only 2.19M SNPs on an older array,
  resulting in 23/41 Tea intake instruments being absent.
- `ebi-a-GCST90018822` (Sakaue 2021) was selected as the best available
  alternative: it covers 24M SNPs and returned all 41 Tea intake instruments,
  but uses ICD-based "chronic renal failure" rather than the eGFR-based CKD
  definition used in the primary analysis. This phenotypic mismatch is the
  primary reason CKD effect estimates are attenuated in the replication.

## Instrument Selection

Genetic instruments were selected using `ieugwasr::tophits()` with the
following parameters (identical to MRAgent defaults):

- P-value threshold: p < 5×10⁻⁸
- LD clumping: R² < 0.001, window 10,000 kb
- Reference panel: 1000 Genomes EUR

Fallback to p < 5×10⁻⁶ was applied when fewer than 11 instruments were
identified at the primary threshold (none required for any exposure in this
analysis).

**Differences from primary analysis:**
The primary analysis additionally applied:
1. PheWAS-based confounder screening (TG, BMI, T2D, SBP; P < 1×10⁻⁵)
2. Exclusion of rs429358 (APOE locus, CAD P = 2.2×10⁻⁹) for Dried fruit,
   Beef, and Lamb/mutton
3. Exclusion of rs1229984 (ADH1B) from all exposures

These filters were not applied in the MRAgent validation to reflect MRAgent's
default (unfiltered) behaviour. Small differences in nSNP between the two
analyses are expected as a result.

## Harmonisation

Exposure and outcome summary statistics were harmonised by:
1. Allele alignment (direct match or strand flip)
2. Removal of palindromic SNPs with EAF within 0.08 of 0.5 (ambiguous strand)
3. Exclusion of SNPs with missing outcome beta

## Statistical Methods

Four MR methods were applied, matching MRAgent's default method set:

| Method | Implementation |
|---|---|
| IVW (random effects) | Multiplicative RE-IVW; SE inflated by √(Cochran's Q / df) when Q/df > 1 |
| Weighted Median | Bowden et al. (2016); SE via bootstrap (n = 1,000) |
| MR-Egger | WLS regression with Bowden et al. (2015) orientation to positive exposure betas |
| Weighted Mode | Hartwig et al. (2017); bandwidth φ = 1, SE via bootstrap (n = 500) |

Multiple testing correction: Benjamini–Hochberg FDR across all 20 exposure–
outcome pairs (10 exposures × 2 outcomes).

## Software

```
R ≥ 4.2
ieugwasr ≥ 0.1.5  (Hemani et al. 2018, doi:10.1101/227452)
dplyr, tidyr (tidyverse)
MRAgent 0.2.5 (Xu et al. 2024, doi:10.1093/bioinformatics/btae619)
```

## References

- Xu W et al. (2024) MRAgent: an automated Mendelian randomisation analysis
  tool powered by large language models. *Bioinformatics* 40(10):btae619.
  https://doi.org/10.1093/bioinformatics/btae619
- Hemani G et al. (2018) The MR-Base platform supports systematic causal
  inference across the human phenome. *eLife* 7:e34408.
- Wuttke M et al. (2019) A catalog of genetic loci associated with kidney
  function from analyses of a million individuals. *Nature Genetics* 51:957–972.
  PMID 31152163
- Sakaue S et al. (2021) A cross-population atlas of genetic associations for
  210 human diseases and traits. *Nature Genetics* 53:1415–1424. PMID 34594039
- Bowden J et al. (2015) Mendelian randomization with invalid instruments.
  *International Journal of Epidemiology* 44(2):512–525.
- Bowden J et al. (2016) Consistent estimation in Mendelian randomization with
  some invalid instruments using a weighted median estimator.
  *Genetic Epidemiology* 40(4):304–314.
- Hartwig FP et al. (2017) Robust inference in summary data Mendelian
  randomization via the zero modal pleiotropy assumption.
  *International Journal of Epidemiology* 46(6):1985–1998.
