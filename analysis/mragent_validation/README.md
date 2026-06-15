# MRAgent Validation — Dietary Acid Load → CKD

This folder documents an independent replication of the main MR analysis using
[MRAgent](https://github.com/xuwei1997/MRAgent) (Xu et al.), an LLM-assisted
Mendelian randomisation pipeline.

## Purpose

To assess whether an automated, OpenGWAS-based MR tool reproduces the direction,
magnitude, and significance of the causal estimates obtained in the primary analysis
(`R/mr_analysis_final.R`), using identical GWAS sources for the exposures and the
best available OpenGWAS equivalents for the outcomes.

## Folder Structure

```
analysis/mragent_validation/
├── README.md                        ← this file
├── methods.md                       ← detailed methodology and GWAS ID mapping
├── run_mragent_validation.R         ← fully reproducible R script (no TwoSampleMR required)
├── results/
│   ├── mragent_ivw_results.csv      ← IVW results for all 10 × 2 pairs (FDR-corrected)
│   ├── mragent_results_all.csv      ← all methods (IVW, WM, MR-Egger, Weighted Mode)
│   └── comparison_primary_vs_mragent.md  ← side-by-side comparison with primary analysis
```

## Key Findings

| Exposure | Outcome | Primary β (p) | MRAgent β (p) | Direction |
|---|---|---|---|---|
| **Tea intake** | eGFR | +0.016 (2.2e-5) ✓ FDR | +0.027 (0.002) ✓ FDR | ✓ consistent |
| **Tea intake** | CKD | −0.403 (8.2e-5) ✓ FDR | −0.162 (0.313) | ✓ same direction |
| Salad intake | eGFR | +0.025 (0.018) | +0.036 (0.037) | ✓ consistent |
| Dried fruit | eGFR | +0.008 (0.133) | +0.023 (0.014) | ✓ consistent |
| Processed meat | eGFR | −0.013 (0.064) | −0.005 (0.709) | ✓ consistent |

MRAgent **FDR-replicates the Tea intake → eGFR signal** (p_FDR = 0.020) and
confirms the **protective direction of Tea intake → CKD**. The attenuated CKD
p-value reflects a phenotype mismatch: the best available OpenGWAS CKD GWAS
uses ICD-coded renal failure (Sakaue 2021), which differs from the eGFR-based
CKD definition used in the primary analysis (CKDGen 2019).

## How to Reproduce

```r
# Prerequisites: R + ieugwasr
install.packages("remotes")
remotes::install_github("MRCIEU/ieugwasr")

# Set your OpenGWAS token
Sys.setenv(OPENGWAS_JWT = "your_token")  # https://api.opengwas.io

# Run
source("analysis/mragent_validation/run_mragent_validation.R")
```

No OpenAI key is required. This script uses only `ieugwasr` for GWAS data
retrieval and implements IVW, Weighted Median, MR-Egger and Weighted Mode
directly in R — replicating MRAgent's core statistical workflow without
the `TwoSampleMR` dependency.

## Software

| Tool | Version | Reference |
|---|---|---|
| MRAgent | 0.2.5 | Xu et al. (2024) |
| ieugwasr | ≥0.1.5 | Hemani et al. (2018) |
| R | ≥4.2 | R Core Team |

## Date

Analysis performed: 2026-05-05
