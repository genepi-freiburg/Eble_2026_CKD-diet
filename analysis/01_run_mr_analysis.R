###############################################################
# Entry Point: Run MR Analysis Pipeline
# 
# Usage:   Rscript analysis/01_run_mr_analysis.R
# 
# Requires:
#   - OPENGWAS_JWT set in environment (see .env.example)
#   - CKDGen GWAS data in data/raw/ (see README for download URLs)
#   - FinnGen R12 data in data/raw/finngen/
#
# Output: all results saved to data/processed/
###############################################################

cat("=== MR Analysis Pipeline ===\n\n")
source("R/mr_analysis_final.R")
