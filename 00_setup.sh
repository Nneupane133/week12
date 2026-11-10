#!/bin/bash
# =============================================================================
# 00_setup.sh
# Install all tools needed for the ARV contamination analysis
# Run this ONCE before anything else
# =============================================================================

set -euo pipefail

echo "========================================"
echo " ARV Contamination Analysis - Setup"
echo "========================================"

# ---------- Check for conda ----------
if ! command -v conda &> /dev/null; then
    echo "[ERROR] conda not found."
    echo "Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
    echo "Then re-run this script."
    exit 1
fi

ENV_NAME="arv_contam"

# ---------- Create environment if it doesn't exist ----------
if conda env list | grep -q "^${ENV_NAME}"; then
    echo "[INFO] Conda environment '${ENV_NAME}' already exists. Skipping creation."
else
    echo "[INFO] Creating conda environment: ${ENV_NAME}"
    conda create -y -n "${ENV_NAME}" python=3.10
fi

# ---------- Install bioinformatics tools ----------
echo "[INFO] Installing bioinformatics tools via conda..."
conda install -y -n "${ENV_NAME}" \
    -c bioconda -c conda-forge \
    sra-tools \
    entrez-direct \
    fastqc \
    fastp \
    bwa-mem2 \
    samtools \
    bedtools

# ---------- Install Python packages ----------
echo "[INFO] Installing Python packages..."
conda run -n "${ENV_NAME}" pip install \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    pysam \
    biopython

# ---------- Verify installs ----------
echo ""
echo "[INFO] Verifying tool installations..."
TOOLS=(prefetch fasterq-dump fastqc fastp bwa-mem2 samtools esearch efetch)
ALL_OK=true
for tool in "${TOOLS[@]}"; do
    if conda run -n "${ENV_NAME}" which "${tool}" &>/dev/null; then
        echo "  [OK] ${tool}"
    else
        echo "  [MISSING] ${tool}"
        ALL_OK=false
    fi
done

echo ""
if [ "$ALL_OK" = true ]; then
    echo "[SUCCESS] All tools installed."
    echo ""
    echo "To activate the environment, run:"
    echo "  conda activate ${ENV_NAME}"
else
    echo "[WARNING] Some tools are missing. Try running setup again or install manually."
fi
