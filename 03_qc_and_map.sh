#!/bin/bash
# =============================================================================
# 03_qc_and_map.sh
# Quality control (fastp) + map all 6 samples to the ARV S1133 reference
# Collect statistics needed for contamination analysis
# =============================================================================

set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
source /apps/x86-64/apps/anaconda_3-2025.12/etc/profile.d/conda.sh
conda activate arv_contam

echo "========================================"
echo " Step 3: QC and Mapping"
echo "========================================"

# ---------- Directories ----------
RAW_DIR="data/raw"
TRIM_DIR="data/trimmed"
BAM_DIR="data/bam"
STATS_DIR="data/stats"
REF="data/reference/ARV_S1133_all_segments.fasta"
THREADS=4

mkdir -p "${TRIM_DIR}" "${BAM_DIR}" "${STATS_DIR}"

SAMPLES=("AL" "S1133" "P1" "PP3" "PPP" "EP")

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}/${SAMPLE}_R1.fastq.gz"
    R2="${RAW_DIR}/${SAMPLE}/${SAMPLE}_R2.fastq.gz"

    # Skip if raw data not present
    if [ ! -f "${R1}" ]; then
        echo "[SKIP] ${SAMPLE}: raw reads not found (${R1})"
        continue
    fi

    echo ""
    echo "-------- Processing ${SAMPLE} --------"

    # ---- A. Adapter trimming with fastp ----
    echo "[INFO] ${SAMPLE}: Trimming adapters with fastp..."
    fastp \
        -i "${R1}" -I "${R2}" \
        -o "${TRIM_DIR}/${SAMPLE}_R1.fastq.gz" \
        -O "${TRIM_DIR}/${SAMPLE}_R2.fastq.gz" \
        --json "${STATS_DIR}/${SAMPLE}_fastp.json" \
        --html "${STATS_DIR}/${SAMPLE}_fastp.html" \
        --thread "${THREADS}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        2>&1 | grep -E "reads|bases|rate|passed"

    # ---- B. Map to ARV S1133 reference ----
    echo "[INFO] ${SAMPLE}: Mapping to S1133 reference..."
    bwa-mem2 mem \
        -t "${THREADS}" \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        "${REF}" \
        "${TRIM_DIR}/${SAMPLE}_R1.fastq.gz" \
        "${TRIM_DIR}/${SAMPLE}_R2.fastq.gz" \
    | samtools sort -@ "${THREADS}" -o "${BAM_DIR}/${SAMPLE}.bam"

    samtools index "${BAM_DIR}/${SAMPLE}.bam"

    # ---- C. Collect mapping statistics ----
    echo "[INFO] ${SAMPLE}: Collecting mapping stats..."

    # flagstat: overall mapping rate
    samtools flagstat "${BAM_DIR}/${SAMPLE}.bam" \
        > "${STATS_DIR}/${SAMPLE}_flagstat.txt"

    # idxstats: per-segment read counts (key for contamination check)
    samtools idxstats "${BAM_DIR}/${SAMPLE}.bam" \
        > "${STATS_DIR}/${SAMPLE}_idxstats.txt"

    # depth: per-position coverage (for correlation analysis)
    # -a = output all positions including 0 depth
    samtools depth -a "${BAM_DIR}/${SAMPLE}.bam" \
        > "${STATS_DIR}/${SAMPLE}_depth.txt"

    # mapq distribution: mapping quality score histogram
    samtools view -F 4 "${BAM_DIR}/${SAMPLE}.bam" \
        | awk '{print $5}' \
        | sort -n | uniq -c \
        > "${STATS_DIR}/${SAMPLE}_mapq_hist.txt"

    echo "[DONE] ${SAMPLE}"
done

echo ""
echo "========================================"
echo " QC and mapping complete."
echo " BAMs:  ${BAM_DIR}/"
echo " Stats: ${STATS_DIR}/"
echo ""
echo " Key files for contamination analysis:"
ls "${STATS_DIR}"/*_flagstat.txt | sort
echo "========================================"
