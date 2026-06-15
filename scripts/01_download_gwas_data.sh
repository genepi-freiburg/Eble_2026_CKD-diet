#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Download outcome GWAS summary statistics required by the MR pipeline.
# Run from repository root:
#     bash scripts/01_download_gwas_data.sh
#
# Total download size: ~1.1 GB
#   - CKDGen eGFR (Wuttke 2019, EA):  ~184 MB
#   - CKDGen CKD  (Wuttke 2019, EA):  ~171 MB
#   - FinnGen R12 N14_CHRONKIDNEYDIS: ~778 MB
# ---------------------------------------------------------------------------
set -euo pipefail

mkdir -p data/raw data/raw/finngen

echo "[1/3] Downloading CKDGen eGFR summary statistics (European ancestry, ~184 MB)..."
wget -c -O data/raw/egfr_EA.txt.gz \
  "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz"

echo "[2/3] Downloading CKDGen CKD summary statistics (European ancestry, ~171 MB)..."
wget -c -O data/raw/ckd_EA.txt.gz \
  "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz"

echo "[3/3] Downloading FinnGen R12 N14_CHRONKIDNEYDIS (~778 MB)..."
# FinnGen R12 manifest (full): https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv
wget -c -O data/raw/finngen/finngen_R12_N14_CHRONKIDNEYDIS.gz \
  "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_N14_CHRONKIDNEYDIS.gz"

echo
echo "Done. Files saved to data/raw/:"
ls -lh data/raw/*.gz data/raw/finngen/*.gz
