#!/usr/bin/env bash
# Run the full pipeline from the repo root. Optional: export GHM_RAW_DIR=/path/for/raw/downloads
set -euo pipefail
cd "$(dirname "$0")"
python -I scripts/01_download_data.py
Rscript scripts/02_extract_may_loads.R
python -I scripts/03_verify_hypoxia_transcription.py
python -I scripts/04_build_dataset_and_regressions.py
python -I scripts/05_bayesian_models.py
python -I scripts/06_figures.py
