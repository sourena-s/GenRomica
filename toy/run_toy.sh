#!/usr/bin/env bash
# Runs the whole pipeline on simulated toy data, then compares it with the
# original code. Usage: bash toy/run_toy.sh   (from any folder)
set -euo pipefail
cd "$(dirname "$0")/.."

Rscript toy/make_toy_data.R
for step in R/01_export_genr_data.R R/02_merge_components.R R/03_clean_components.R \
            R/04_compute_scores.R R/05_run_regressions.R; do
  echo "=== $step"
  Rscript "$step" config/toy.R
done
echo "=== Comparing with the original code"
Rscript toy/compare_with_original.R
