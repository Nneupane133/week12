# ARV Contamination Analysis
**BioProject:** PRJNA1217582  
**Paper:** Avian Diseases, Vol. 70, No. 1, p. 23 (2026)  
**Question:** Are the S1133-like reads in the Alabama (AL) field isolate of Avian Reovirus genuine co-infection — or index-hopping contamination from the sequenced S1133 control sample?

---

## Background

The original paper sequenced an AL field isolate and found it contained two populations: one similar to the reference vaccine strain (S1133) and one divergent field strain. They concluded the mixture persisted even after 3 rounds of plaque purification.

**This analysis re-tests that conclusion.** When multiple samples are sequenced on the same Illumina run, index-hopping (barcode switching, ~0.1–1% rate) can transfer reads between samples. If the AL isolate and the S1133 control were on the same MiSeq run, the "S1133-like" reads in AL could simply be contamination — which would mean the isolate was never truly mixed, and the entire paper's conclusion rests on a technical artifact.

---

## Samples

| Label  | SRX Accession  | Description                             |
|--------|----------------|-----------------------------------------|
| AL     | SRX27528194    | AL field isolate ← **key sample**       |
| S1133  | SRX30741467    | Pure S1133 control ← **potential contaminant source** |
| P1     | SRX30741463    | Mixed AL-S1133, 1st CELi passage        |
| PP3    | SRX30741464    | Mixed AL-S1133, 3× plaque purification  |
| PPP    | SRX30741465    | Mixed AL-S1133, post-purification CELi  |
| EP     | SRX30741466    | Mixed AL-S1133, embryo passage in CELi  |

---

## Four Statistical Tests

### Test 1 — Mapping fraction comparison
Map all reads to the ARV S1133 reference. If the AL isolate has <1% of its reads mapping to S1133 (vs. ~95%+ for the pure S1133 control), that level is consistent with contamination and not co-infection.

### Test 2 — Index-hopping rate estimation  
Model what fraction of AL reads would be S1133-derived under different hopping rates (0.1%, 0.5%, 1%). If the observed S1133-like fraction in AL matches any of these modeled rates, contamination is a plausible and sufficient explanation.

### Test 3 — Per-position Pearson correlation  
If contamination: the coverage depth across the S1133 genome in the AL isolate should closely mirror the pure S1133 control (same library prep biases → same coverage bumps and dips → high r).  
If true co-infection: the two coverage profiles should be independent.

### Test 4 — Gini coefficient (coverage uniformity)  
Contamination tends to accumulate in highly conserved / high-identity regions, producing uneven coverage (high Gini). True viral replication drives more uniform genome coverage. If AL's Gini is substantially higher than S1133's Gini, that supports contamination.

---

## Expected Outputs

```
results/
├── contamination_report.txt          # Text summary of all 4 tests with verdict
└── figures/
    └── contamination_analysis.png    # 4-panel figure (mapping %, model, correlation, Gini)
```

---

## How to Run

```bash
# 1. Install tools (once)
bash 00_setup.sh

# 2. Run full pipeline
bash run_all.sh

# Or step by step:
bash 01_download_sra.sh       # ~5 min (data is small, ~25 MB total)
bash 02_get_reference.sh      # ~2 min
bash 03_qc_and_map.sh         # ~5 min
python3 04_contamination_analysis.py  # ~1 min
```

---

## Interpretation Guide

| Outcome | Contamination Supported | Co-infection Supported |
|---------|------------------------|------------------------|
| Test 1 (mapping %) | AL < 1% | AL > 5% |
| Test 2 (hop model) | Observed ≈ predicted at 0.1–1% | Observed >> all hop rate predictions |
| Test 3 (Pearson r) | r > 0.85 | r < 0.50 |
| Test 4 (Gini) | AL Gini >> S1133 Gini | AL Gini ≈ S1133 Gini |

**If contamination is supported:** The original paper's conclusion is an artifact of the sequencing run. The isolate may be a pure divergent field strain, not a mixed population. This has major implications — it would mean plaque purification worked, and the S1133-like contigs were sequencing noise.

**If co-infection is supported:** The original paper's conclusion stands, and the question becomes *why* plaque purification failed — a separate biological question about growth advantage of mixed infections.

---

## Tools Required
- SRA Toolkit (prefetch, fasterq-dump)
- Entrez Direct (esearch, efetch)
- fastp
- BWA-MEM2
- samtools
- Python 3.10+: pandas, numpy, matplotlib, seaborn, scipy, pysam, biopython

All installed automatically via `00_setup.sh` using conda.
