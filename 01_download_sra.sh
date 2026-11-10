#!/bin/bash
# =============================================================================
# 01_download_sra.sh
# Download all 6 SRA runs from BioProject PRJNA1217582
#
# Sample map:
#   AL        = SRX27528194  → AL field isolate (THE KEY SAMPLE)
#   S1133     = SRX30741467  → Pure S1133 control strain (potential contaminator)
#   P1        = SRX30741463  → Mixed AL-S1133, 1st CELi passage
#   PP3       = SRX30741464  → Mixed AL-S1133, after 3 rounds plaque purification
#   PPP       = SRX30741465  → Mixed AL-S1133, post-purification CELi passage
#   EP        = SRX30741466  → Mixed AL-S1133, embryo passage in CELi
# =============================================================================

set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
source /apps/x86-64/apps/anaconda_3-2025.12/etc/profile.d/conda.sh
conda activate arv_contam

echo "========================================"
echo " Step 1: Download SRA Data"
echo "========================================"

# ---------- Directories ----------
RAW_DIR="data/raw"
mkdir -p "${RAW_DIR}"

# ---------- Sample definitions ----------
declare -A SAMPLES
SAMPLES["AL"]="SRX27528194"
SAMPLES["S1133"]="SRX30741467"
SAMPLES["P1"]="SRX30741463"
SAMPLES["PP3"]="SRX30741464"
SAMPLES["PPP"]="SRX30741465"
SAMPLES["EP"]="SRX30741466"

# ---------- Download each sample ----------
for SAMPLE in "${!SAMPLES[@]}"; do
    SRX="${SAMPLES[$SAMPLE]}"
    OUT_DIR="${RAW_DIR}/${SAMPLE}"
    mkdir -p "${OUT_DIR}"

    echo ""
    echo "[INFO] Downloading ${SAMPLE} (${SRX})..."

    # prefetch downloads the .sra file
    prefetch "${SRX}" --output-directory "${OUT_DIR}" --progress

    # fasterq-dump converts .sra to paired FASTQ
    # prefetch saves files by run accession (SRR*/SRR*.sra), not by SRX,
    # so convert every .sra file it downloaded for this experiment
    # --split-files = separate R1 and R2
    shopt -s nullglob
    SRA_FILES=("${OUT_DIR}"/SRR*/*.sra)
    if [ ${#SRA_FILES[@]} -eq 0 ]; then
        echo "[ERROR] No .sra files found in ${OUT_DIR} after prefetch" >&2
        exit 1
    fi
    for SRA in "${SRA_FILES[@]}"; do
        fasterq-dump "${SRA}" \
            --outdir "${OUT_DIR}" \
            --split-files \
            --threads 4 \
            --progress
    done

    # Combine all runs into consistent names: SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz
    # (concatenating handles experiments with more than one SRR run)
    echo "[INFO] Merging and compressing ${SAMPLE} FASTQs..."
    R1_FILES=("${OUT_DIR}"/SRR*_1.fastq)
    R2_FILES=("${OUT_DIR}"/SRR*_2.fastq)
    if [ ${#R1_FILES[@]} -eq 0 ]; then
        echo "[ERROR] fasterq-dump produced no *_1.fastq files for ${SAMPLE}" >&2
        exit 1
    fi
    cat "${R1_FILES[@]}" | gzip > "${OUT_DIR}/${SAMPLE}_R1.fastq.gz"
    rm -f "${R1_FILES[@]}"
    if [ ${#R2_FILES[@]} -gt 0 ]; then
        cat "${R2_FILES[@]}" | gzip > "${OUT_DIR}/${SAMPLE}_R2.fastq.gz"
        rm -f "${R2_FILES[@]}"
    else
        echo "[WARN] ${SAMPLE} has no R2 reads (single-end data?)"
    fi
    shopt -u nullglob

    echo "[DONE] ${SAMPLE}: ${OUT_DIR}/${SAMPLE}_R1.fastq.gz"
done

echo ""
echo "========================================"
echo " Download complete. Files in: ${RAW_DIR}/"
echo " Directory structure:"
find "${RAW_DIR}" -name "*.fastq.gz" | sort
echo "========================================"
