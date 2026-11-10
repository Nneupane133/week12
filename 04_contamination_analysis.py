#!/usr/bin/env python3
"""
04_contamination_analysis.py
============================
Tests whether the S1133-like reads in the AL isolate are contamination
from the pure S1133 control sample (index-hopping) or true co-infection.

FOUR TESTS:
  Test 1 – Mapping fraction comparison
           How many reads in each sample map to the S1133 reference?
           Contamination expected: AL isolate has very low % (<1%)
           Co-infection expected:  AL isolate has meaningful % (>5%)

  Test 2 – Index-hopping rate estimation
           If sequenced on the same MiSeq run, expected contamination from
           index-hopping is ~0.1–1% of the total library.
           We calculate what % of AL reads would be S1133 under that model.

  Test 3 – Per-position coverage Pearson correlation
           Contamination: AL's S1133-like reads should mirror the S1133
                          control's coverage pattern (same biases, same library)
           Co-infection:  AL's coverage should be independent of control

  Test 4 – Coverage uniformity (Gini coefficient)
           Contamination: uneven coverage, peaking at conserved regions
           Co-infection:  uniform coverage across all 10 segments

Outputs:
  results/contamination_report.txt  — text summary
  results/figures/                  — publication-ready figures
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr

# ---------- Paths ----------
STATS_DIR  = "data/stats"
OUT_DIR    = "results"
FIG_DIR    = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLES = {
    "AL":    "AL Isolate (field)",
    "S1133": "S1133 Control (pure)",
    "P1":    "1st CELi Passage",
    "PP3":   "3× Plaque Purification",
    "PPP":   "Post-Purif. CELi",
    "EP":    "Embryo Passage",
}

# Index-hopping rates to model (Illumina MiSeq typical range)
HOP_RATES = [0.001, 0.005, 0.01]  # 0.1%, 0.5%, 1%

# ============================================================
# Helper functions
# ============================================================

def gini(array):
    """Gini coefficient of a 1D array (measure of inequality/unevenness)."""
    array = np.sort(np.array(array, dtype=float))
    n = len(array)
    if n == 0 or array.sum() == 0:
        return np.nan
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * array) / (n * array.sum())) - (n + 1) / n


def read_flagstat(sample):
    """Parse samtools flagstat output → dict with key stats."""
    path = os.path.join(STATS_DIR, f"{sample}_flagstat.txt")
    if not os.path.exists(path):
        return None
    result = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if "in total" in line:
                result["total"] = int(parts[0])
            elif "mapped (" in line and "mate" not in line:
                result["mapped"] = int(parts[0])
    if "total" in result and "mapped" in result:
        result["pct_mapped"] = 100 * result["mapped"] / result["total"] if result["total"] else 0
    return result


def read_depth(sample):
    """Read samtools depth file → dict keyed by segment name."""
    path = os.path.join(STATS_DIR, f"{sample}_depth.txt")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["segment", "pos", "depth"])
    return df


def read_mapq(sample):
    """Read mapping quality histogram."""
    path = os.path.join(STATS_DIR, f"{sample}_mapq_hist.txt")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                rows.append({"count": int(parts[0]), "mapq": int(parts[1])})
    return pd.DataFrame(rows)


# ============================================================
# TEST 1 — Mapping fraction
# ============================================================

def test1_mapping_fraction(report):
    report.append("\n" + "="*60)
    report.append("TEST 1: Mapping Fraction to S1133 Reference")
    report.append("="*60)
    report.append(f"{'Sample':<12} {'Total Reads':>14} {'Mapped to S1133':>16} {'% Mapped':>10}")
    report.append("-"*55)

    data = {}
    for s in SAMPLES:
        fs = read_flagstat(s)
        if fs is None:
            report.append(f"  {s:<12} [data not found - run 03_qc_and_map.sh first]")
            continue
        data[s] = fs
        report.append(f"  {s:<12} {fs['total']:>14,} {fs['mapped']:>16,} {fs['pct_mapped']:>9.2f}%")

    # Interpretation
    report.append("")
    if "AL" in data and "S1133" in data:
        al_pct  = data["AL"]["pct_mapped"]
        s1_pct  = data["S1133"]["pct_mapped"]
        ratio   = al_pct / s1_pct if s1_pct > 0 else 0

        report.append(f"AL isolate maps {al_pct:.2f}% of reads to S1133 reference.")
        report.append(f"Pure S1133 control maps {s1_pct:.2f}% of reads to S1133 reference.")
        report.append(f"Ratio (AL / S1133): {ratio:.4f}")
        report.append("")
        if ratio < 0.01:
            report.append("→ CONTAMINATION LIKELY: AL maps <1% of what S1133 control maps.")
            report.append("  This is consistent with index-hopping contamination.")
        elif ratio < 0.05:
            report.append("→ AMBIGUOUS: AL maps 1-5% of S1133 control. Could be low-level")
            report.append("  co-infection or contamination. See Tests 3 and 4.")
        else:
            report.append("→ CO-INFECTION LIKELY: AL maps >5% relative to S1133 control.")
            report.append("  This exceeds typical index-hopping rates (0.1-1%).")
    return data


# ============================================================
# TEST 2 — Index-hopping rate estimation
# ============================================================

def test2_index_hopping(flagstat_data, report):
    report.append("\n" + "="*60)
    report.append("TEST 2: Index-Hopping Rate Estimation")
    report.append("="*60)
    report.append("If the S1133 control and AL isolate were sequenced on the")
    report.append("same MiSeq run, index-hopping can transfer ~0.1-1% of reads")
    report.append("from one sample's index to another's.")
    report.append("")

    if "AL" not in flagstat_data or "S1133" not in flagstat_data:
        report.append("  [Missing data for AL or S1133]")
        return

    s1133_total = flagstat_data["S1133"]["total"]
    al_total    = flagstat_data["AL"]["total"]
    al_mapped   = flagstat_data["AL"]["mapped"]
    al_obs_pct  = flagstat_data["AL"]["pct_mapped"]

    report.append(f"S1133 control total reads: {s1133_total:,}")
    report.append(f"AL isolate total reads:    {al_total:,}")
    report.append(f"AL reads mapped to S1133:  {al_mapped:,} ({al_obs_pct:.3f}%)")
    report.append("")
    report.append(f"{'Hop Rate':>12} {'Expected S1133 reads in AL':>28} {'Expected %':>12} {'Observed %':>12} {'Verdict':>14}")
    report.append("-"*80)

    for rate in HOP_RATES:
        expected_reads = s1133_total * rate
        expected_pct   = 100 * expected_reads / al_total if al_total > 0 else 0
        verdict = "MATCHES" if abs(al_obs_pct - expected_pct) < 0.5 else \
                  ("BELOW" if al_obs_pct < expected_pct else "EXCEEDS")
        report.append(
            f"  {rate*100:.1f}%        {expected_reads:>28,.0f} {expected_pct:>11.3f}% "
            f"{al_obs_pct:>11.3f}%  {verdict:>14}"
        )

    report.append("")
    report.append("If 'MATCHES' at any hop rate → contamination plausible.")
    report.append("If 'EXCEEDS' all hop rates   → true co-infection more likely.")


# ============================================================
# TEST 3 — Per-position coverage Pearson correlation
# ============================================================

def test3_coverage_correlation(report):
    report.append("\n" + "="*60)
    report.append("TEST 3: Per-Position Coverage Pearson Correlation")
    report.append("="*60)
    report.append("Contamination prediction: AL's coverage pattern of S1133-like")
    report.append("reads should closely mirror the pure S1133 control (same library")
    report.append("preparation biases → high positive correlation, r ≈ 0.9+).")
    report.append("Co-infection prediction: independent coverage, lower correlation.")
    report.append("")

    al_depth  = read_depth("AL")
    s1_depth  = read_depth("S1133")

    if al_depth is None or s1_depth is None:
        report.append("  [Depth files not found - run 03_qc_and_map.sh first]")
        return None, None

    # Merge on segment and position
    merged = pd.merge(al_depth, s1_depth,
                      on=["segment", "pos"],
                      suffixes=("_AL", "_S1133"))

    if len(merged) < 100:
        report.append("  [Insufficient overlapping positions for correlation]")
        return None, None

    r, p = pearsonr(merged["depth_AL"], merged["depth_S1133"])
    report.append(f"Pearson r = {r:.4f}  (p = {p:.2e})")
    report.append(f"Positions compared: {len(merged):,}")
    report.append("")
    if r > 0.85:
        report.append("→ HIGH CORRELATION (r > 0.85): Coverage patterns are very similar.")
        report.append("  CONTAMINATION SUPPORTED — AL's S1133-like reads mirror control.")
    elif r > 0.5:
        report.append("→ MODERATE CORRELATION: Ambiguous. Consider other tests.")
    else:
        report.append("→ LOW CORRELATION (r < 0.5): Coverage patterns are independent.")
        report.append("  CO-INFECTION SUPPORTED — AL replicates S1133 on its own trajectory.")

    return merged, r


# ============================================================
# TEST 4 — Gini coefficient (coverage uniformity)
# ============================================================

def test4_gini(report):
    report.append("\n" + "="*60)
    report.append("TEST 4: Coverage Uniformity (Gini Coefficient)")
    report.append("="*60)
    report.append("Contamination: S1133-like reads in AL are biased toward")
    report.append("conserved/high-identity regions → uneven coverage, high Gini.")
    report.append("Co-infection:  replication-driven coverage is relatively uniform,")
    report.append("               Gini similar to what pure S1133 control shows.")
    report.append("")
    report.append(f"{'Sample':<12} {'Mean Depth':>12} {'Gini':>8}  Interpretation")
    report.append("-"*60)

    gini_vals = {}
    for s in SAMPLES:
        df = read_depth(s)
        if df is None:
            report.append(f"  {s:<12} [not found]")
            continue
        depths = df["depth"].values
        g = gini(depths)
        mean_d = depths.mean()
        gini_vals[s] = g
        interp = "uniform" if g < 0.4 else ("moderate" if g < 0.7 else "uneven")
        report.append(f"  {s:<12} {mean_d:>12.1f} {g:>8.4f}  {interp}")

    report.append("")
    if "AL" in gini_vals and "S1133" in gini_vals:
        al_g  = gini_vals["AL"]
        s1_g  = gini_vals["S1133"]
        diff  = al_g - s1_g
        report.append(f"AL Gini = {al_g:.4f}, S1133 Gini = {s1_g:.4f}, Difference = {diff:.4f}")
        if diff > 0.15:
            report.append("→ AL coverage is MORE UNEVEN than S1133 control.")
            report.append("  CONTAMINATION SUPPORTED — reads accumulate at conserved regions.")
        elif abs(diff) < 0.05:
            report.append("→ Similar uniformity to S1133 control.")
            report.append("  CO-INFECTION SUPPORTED — independent, uniform replication.")
        else:
            report.append("→ Moderate difference. Ambiguous.")
    return gini_vals


# ============================================================
# FIGURES
# ============================================================

def make_figures(flagstat_data, depth_merged, pearson_r, gini_vals):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    colors = {
        "AL":    "#e63946",   # red — the key sample
        "S1133": "#457b9d",   # blue — control
        "P1":    "#2a9d8f",
        "PP3":   "#e9c46a",
        "PPP":   "#f4a261",
        "EP":    "#a8dadc",
    }

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(
        "ARV AL Isolate — Contamination vs. Co-infection Analysis\n"
        "(S1133-like reads in AL isolate: index-hopping or true mixed infection?)",
        fontsize=13, fontweight="bold", y=0.98
    )
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

    # ---- Panel A: Mapping fraction bar chart ----
    ax_a = fig.add_subplot(gs[0, 0])
    if flagstat_data:
        samps  = [s for s in SAMPLES if s in flagstat_data]
        pcts   = [flagstat_data[s]["pct_mapped"] for s in samps]
        bars   = ax_a.bar(samps, pcts,
                          color=[colors.get(s, "gray") for s in samps],
                          edgecolor="white", linewidth=0.5)
        ax_a.axhline(1.0, ls="--", color="gray", lw=1,
                     label="1% threshold (max index-hopping)")
        ax_a.set_ylabel("% Reads Mapped to S1133 Reference")
        ax_a.set_title("A  Mapping Fraction per Sample")
        ax_a.legend(fontsize=8)
        for bar, pct in zip(bars, pcts):
            ax_a.text(bar.get_x() + bar.get_width()/2,
                      bar.get_height() + 0.05,
                      f"{pct:.2f}%", ha="center", va="bottom", fontsize=8)

    # ---- Panel B: Index-hopping model ----
    ax_b = fig.add_subplot(gs[0, 1])
    if flagstat_data and "AL" in flagstat_data and "S1133" in flagstat_data:
        s1_total = flagstat_data["S1133"]["total"]
        al_total = flagstat_data["AL"]["total"]
        al_obs   = flagstat_data["AL"]["pct_mapped"]

        rates    = np.linspace(0.0001, 0.02, 200)
        exp_pcts = [100 * s1_total * r / al_total for r in rates]

        ax_b.plot(rates * 100, exp_pcts, "k-", lw=1.5, label="Expected (model)")
        ax_b.axhline(al_obs, color=colors["AL"], lw=2, ls="--",
                     label=f"Observed in AL ({al_obs:.3f}%)")
        ax_b.fill_between(rates * 100, exp_pcts, al_obs,
                          where=[e < al_obs for e in exp_pcts],
                          alpha=0.15, color=colors["AL"],
                          label="Excess above model")
        ax_b.set_xlabel("Assumed Index-Hopping Rate (%)")
        ax_b.set_ylabel("Expected S1133 Reads in AL (%)")
        ax_b.set_title("B  Index-Hopping Rate Model")
        ax_b.legend(fontsize=8)

    # ---- Panel C: Coverage correlation scatter ----
    ax_c = fig.add_subplot(gs[1, 0])
    if depth_merged is not None and len(depth_merged) > 0:
        ax_c.scatter(depth_merged["depth_S1133"],
                     depth_merged["depth_AL"],
                     alpha=0.2, s=2, color="#264653", rasterized=True)
        if pearson_r is not None:
            ax_c.set_title(f"C  Coverage Correlation  (r = {pearson_r:.3f})")
        else:
            ax_c.set_title("C  Coverage Correlation")
        ax_c.set_xlabel("S1133 Control — Depth per Position")
        ax_c.set_ylabel("AL Isolate — Depth per Position")
        # Fit line
        m, b, _, _, _ = stats.linregress(depth_merged["depth_S1133"],
                                          depth_merged["depth_AL"])
        xmax = depth_merged["depth_S1133"].max()
        ax_c.plot([0, xmax], [b, m * xmax + b], "r-", lw=1.5, alpha=0.7)
        note = "High r → contamination" if (pearson_r or 0) > 0.7 else "Low r → co-infection"
        ax_c.text(0.05, 0.92, note, transform=ax_c.transAxes,
                  fontsize=8, color="red" if (pearson_r or 0) > 0.7 else "green")

    # ---- Panel D: Gini coefficients ----
    ax_d = fig.add_subplot(gs[1, 1])
    if gini_vals:
        samps = [s for s in SAMPLES if s in gini_vals]
        gs_v  = [gini_vals[s] for s in samps]
        ax_d.bar(samps, gs_v,
                 color=[colors.get(s, "gray") for s in samps],
                 edgecolor="white", linewidth=0.5)
        ax_d.axhline(0.4, ls="--", color="gray", lw=1,
                     label="Uniformity threshold (0.4)")
        ax_d.set_ylabel("Gini Coefficient (0 = uniform, 1 = uneven)")
        ax_d.set_title("D  Coverage Uniformity (Gini)")
        ax_d.legend(fontsize=8)
        for bar, g in zip(ax_d.patches, gs_v):
            ax_d.text(bar.get_x() + bar.get_width()/2,
                      bar.get_height() + 0.005,
                      f"{g:.3f}", ha="center", va="bottom", fontsize=8)

    # Save
    fig_path = os.path.join(FIG_DIR, "contamination_analysis.png")
    plt.savefig(fig_path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[SAVED] Figure: {fig_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    report = []
    report.append("=" * 60)
    report.append("ARV AL ISOLATE — CONTAMINATION ANALYSIS REPORT")
    report.append("BioProject: PRJNA1217582")
    report.append("Hypothesis: S1133-like reads in AL isolate are index-hopping")
    report.append("            contamination, NOT true co-infection.")
    report.append("=" * 60)

    flagstat_data = test1_mapping_fraction(report)
    test2_index_hopping(flagstat_data, report)
    depth_merged, pearson_r = test3_coverage_correlation(report)
    gini_vals = test4_gini(report)

    # --- Overall verdict ---
    report.append("\n" + "=" * 60)
    report.append("OVERALL INTERPRETATION")
    report.append("=" * 60)
    report.append("")
    report.append("Contamination is SUPPORTED if:")
    report.append("  ✓ AL maps <1% of reads to S1133")
    report.append("  ✓ AL mapping % matches index-hopping model (<1% hop rate)")
    report.append("  ✓ Per-position coverage correlation between AL and S1133 is high (r>0.85)")
    report.append("  ✓ AL Gini is higher than S1133 Gini (uneven coverage)")
    report.append("")
    report.append("Co-infection is SUPPORTED if:")
    report.append("  ✓ AL maps 5-50% of reads to S1133")
    report.append("  ✓ AL mapping % exceeds all index-hopping models")
    report.append("  ✓ Per-position coverage correlation is low (r<0.5)")
    report.append("  ✓ AL Gini is similar to S1133 Gini (uniform coverage)")

    # Print and save report
    report_text = "\n".join(report)
    print(report_text)

    report_path = os.path.join(OUT_DIR, "contamination_report.txt")
    with open(report_path, "w") as f:
        f.write(report_text)
    print(f"\n[SAVED] Report: {report_path}")

    # Make figures
    make_figures(flagstat_data, depth_merged, pearson_r, gini_vals)

    print("\n[DONE] Analysis complete.")
    print(f"  Report: {report_path}")
    print(f"  Figure: {FIG_DIR}/contamination_analysis.png")


if __name__ == "__main__":
    main()
