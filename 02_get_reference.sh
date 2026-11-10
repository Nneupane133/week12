#!/bin/bash
# =============================================================================
# 02_get_reference.sh
# Download the ARV S1133 reference genome (all 10 dsRNA segments) from NCBI
# and prepare it for BWA-MEM2 mapping
#
# ARV S1133 is the classic vaccine/reference strain of Avian Orthoreovirus.
# Its genome has 10 double-stranded RNA segments: L1-L3, M1-M3, S1-S4.
# We use this as the reference to test whether AL isolate reads map to it.
#
# Reference: GenBank KF741756-KF741765, one complete sequence per segment
# ("Avian orthoreovirus strain S1133 segment X, complete sequence")
# =============================================================================

set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
source /apps/x86-64/apps/anaconda_3-2025.12/etc/profile.d/conda.sh
conda activate arv_contam

echo "========================================"
echo " Step 2: Download ARV S1133 Reference"
echo "========================================"

REF_DIR="data/reference"
mkdir -p "${REF_DIR}"

REF_FA="${REF_DIR}/ARV_S1133_all_segments.fasta"

# ---------- S1133 segment accessions (GenBank KF741756-KF741765) ----------
ACCESSIONS=(
    "KF741756"   # L1  3959 bp
    "KF741757"   # L2  3830 bp
    "KF741758"   # L3  3907 bp
    "KF741759"   # M1  2283 bp
    "KF741760"   # M2  2158 bp
    "KF741761"   # M3  1996 bp
    "KF741762"   # S1  1643 bp (sigma C gene - KEY antigen)
    "KF741763"   # S2  1324 bp
    "KF741764"   # S3  1202 bp
    "KF741765"   # S4  1192 bp
)
EXPECTED=${#ACCESSIONS[@]}

# ---------- Download all segments in one request ----------
echo "[INFO] Fetching ${EXPECTED} S1133 segments from NCBI..."
ID_LIST=$(IFS=,; echo "${ACCESSIONS[*]}")

# NCBI occasionally drops connections, so retry a few times
for ATTEMPT in 1 2 3; do
    if efetch -db nuccore -id "${ID_LIST}" -format fasta > "${REF_FA}.tmp" \
        && [ "$(grep -c '^>' "${REF_FA}.tmp" || true)" -eq "${EXPECTED}" ]; then
        break
    fi
    echo "[WARN] Attempt ${ATTEMPT} failed or incomplete; retrying in 10s..."
    sleep 10
done

# ---------- Verify we got exactly the 10 expected S1133 segments ----------
SEQ_COUNT=$(grep -c "^>" "${REF_FA}.tmp" || true)
if [ "${SEQ_COUNT}" -ne "${EXPECTED}" ]; then
    echo "[ERROR] Expected ${EXPECTED} sequences but got ${SEQ_COUNT}." >&2
    rm -f "${REF_FA}.tmp"
    exit 1
fi
for ACC in "${ACCESSIONS[@]}"; do
    if ! grep -q "^>${ACC}" "${REF_FA}.tmp"; then
        echo "[ERROR] ${ACC} missing from download." >&2
        rm -f "${REF_FA}.tmp"
        exit 1
    fi
done
if grep "^>" "${REF_FA}.tmp" | grep -qv "S1133"; then
    echo "[ERROR] A downloaded record is not labeled S1133:" >&2
    grep "^>" "${REF_FA}.tmp" | grep -v "S1133" >&2
    rm -f "${REF_FA}.tmp"
    exit 1
fi
mv "${REF_FA}.tmp" "${REF_FA}"

# ---------- Show what we have ----------
echo ""
echo "[INFO] Reference sequences (${SEQ_COUNT}):"
grep "^>" "${REF_FA}"

# ---------- Index reference for BWA-MEM2 ----------
echo ""
echo "[INFO] Indexing reference with BWA-MEM2..."
bwa-mem2 index "${REF_FA}"

# ---------- Index for samtools ----------
echo "[INFO] Creating samtools fai index..."
samtools faidx "${REF_FA}"

echo ""
echo "[INFO] Segment lengths:"
cut -f1,2 "${REF_FA}.fai"

echo ""
echo "========================================"
echo " Reference ready: ${REF_FA}"
echo "========================================"
