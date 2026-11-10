#!/bin/bash
# =============================================================================
# 04_run_analysis.sh
# PBS wrapper: run the contamination analysis (04_contamination_analysis.py)
# Reads:   data/stats/*  (from 03_qc_and_map.sh)
# Writes:  results/contamination_report.txt, results/figures/
# =============================================================================

set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
source /apps/x86-64/apps/anaconda_3-2025.12/etc/profile.d/conda.sh
conda activate arv_contam

echo "========================================"
echo " Step 4: Contamination Analysis"
echo "========================================"

python 04_contamination_analysis.py

echo ""
echo "[INFO] Outputs:"
ls -l results/ results/figures/
