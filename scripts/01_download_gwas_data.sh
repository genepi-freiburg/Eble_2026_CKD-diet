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
#
# Works on Linux and macOS: uses wget if available, otherwise falls back to
# curl. Both support resuming a partially downloaded file, so re-running the
# script after an interrupted download continues where it left off.
# ---------------------------------------------------------------------------
set -euo pipefail

# --- Select a download tool (wget preferred, curl as fallback) -------------
if command -v wget >/dev/null 2>&1; then
  DL_TOOL="wget"
elif command -v curl >/dev/null 2>&1; then
  DL_TOOL="curl"
else
  echo "ERROR: neither 'wget' nor 'curl' is installed." >&2
  echo "Install one of them and re-run, e.g.:" >&2
  echo "  macOS (Homebrew):  brew install wget" >&2
  echo "  Debian/Ubuntu:     sudo apt-get install wget" >&2
  exit 1
fi
echo "Using download tool: ${DL_TOOL}"

# download <url> <output_path>  — resumable, follows redirects
download() {
  url="$1"; out="$2"
  if [ "${DL_TOOL}" = "wget" ]; then
    wget -c -O "${out}" "${url}"
  else
    # -L follow redirects, -C - resume, --fail error on HTTP >= 400
    curl -L -C - --fail -o "${out}" "${url}"
  fi
}

mkdir -p data/raw data/raw/finngen

echo "[1/3] Downloading CKDGen eGFR summary statistics (European ancestry, ~184 MB)..."
download \
  "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz" \
  "data/raw/egfr_EA.txt.gz"

echo "[2/3] Downloading CKDGen CKD summary statistics (European ancestry, ~171 MB)..."
download \
  "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz" \
  "data/raw/ckd_EA.txt.gz"

echo "[3/3] Downloading FinnGen R12 N14_CHRONKIDNEYDIS (~778 MB)..."
# FinnGen R12 manifest (full): https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv
download \
  "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_N14_CHRONKIDNEYDIS.gz" \
  "data/raw/finngen/finngen_R12_N14_CHRONKIDNEYDIS.gz"

echo
echo "Done. Files saved to data/raw/:"
ls -lh data/raw/*.gz data/raw/finngen/*.gz
