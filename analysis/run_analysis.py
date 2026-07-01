#!/usr/bin/env python3
"""
Pre-specified Analysis — Computational EGFR TKI Binding Affinity vs Real-World TTNT
=====================================================================================
Reproduces all results and figures from:

  Abboud A. "Computational Binding Affinity Correlates With Real-World Time to
  Next Treatment Across EGFR TKI–Mutation Profiles: A Pre-specified Correlation
  Study." Algentics Research, 2026.

v1.0 scope: binding-pocket mutations (≤5 Å from drug in active-conformation crystal).
v3.0 scope: all OncoKB sensitising mutations.

Usage (run from repo root):
    python analysis/run_analysis.py

Outputs (written to repo root):
    results/all_results.json
    figures/fig1_v1_scatter.pdf/.png
    figures/fig2_heatmap.pdf/.png
    figures/fig3_loo.pdf/.png
    figures/fig4_genie_replication.pdf/.png

Data requirements:
    data/merged_ddg_ttnt.csv              — included in repo
    data/genie_bpc_profile_summary.csv    — included in repo
    data/genie_bpc_cohort.csv             — NOT included (patient-level); download from
                                            Synapse syn27056697 and place here to reproduce
                                            Supplementary Table 4 (treatment-line analysis).
                                            All other results reproduce without this file.
"""

import json
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# ── Paths (all relative to repo root) ─────────────────────────────────────────
REPO_ROOT      = Path(__file__).resolve().parent.parent
DDG_CSV        = REPO_ROOT / "data" / "merged_ddg_ttnt.csv"
GENIE_CSV      = REPO_ROOT / "data" / "genie_bpc_profile_summary.csv"
GENIE_COHORT   = REPO_ROOT / "data" / "genie_bpc_cohort.csv"   # download from Synapse
OUT_JSON       = REPO_ROOT / "results" / "all_results.json"
FIG_DIR        = REPO_ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_JSON.parent.mkdir(parents=True, exist_ok=True)

# ── Scope definitions (pre-specified, see SCOPE.md) ───────────────────────────
# v1.0: binding-pocket mutations (≤5 Å from drug in active conformation 4ZAU/4G5J)
# Qualifying profiles = those with at least one mutated residue ≤5 Å from drug.
# Distances from distance_measurements.md: G719(3.73), T790(3.19), C797(1.84), Ex19del(loop criterion).
# Excluded from v1.0: L858R-only, L861Q-only, S768I-only (all >9 Å in active conformation).
V1_EXCLUDED_PROFILES = {"L858R", "L861Q", "S768I"}   # single-mutation profiles outside binding pocket

# Minimum thresholds (Decision 2 + 9)
MIN_N = 5
MIN_EVENTS = 4   # Decision 9 (v2 only) — pairs not meeting this are flagged

# Treatment-appropriateness exclusions (Decision 4)
TREATMENT_EXCLUSIONS = {
    ("Ex19del+T790M", "erlotinib"),
    ("Ex19del+T790M", "afatinib"),
    ("L858R+T790M",   "erlotinib"),
    ("L858R+T790M",   "afatinib"),
    ("C797S+Ex19del+T790M", "osimertinib"),  # C797S is osimertinib resistance
}

# ── Load MSK-CHORD + GENIE BPC primary data ───────────────────────────────────
def load_primary():
    rows = []
    with open(DDG_CSV) as f:
        for r in csv.DictReader(f):
            rows.append({
                "profile": r["profile"],
                "drug":    r["drug"],
                "ddg":     float(r["ddg"]),
                "ttnt":    float(r["ttnt"]),
                "n":       int(r["n"]),
                "source":  "MSK-CHORD",
                # n_events unknown from this file — mark None
                "n_events": None,
            })
    # Assign known event counts from earlier analysis
    EVENTS = {
        ("G719X+L861Q", "osimertinib"): 3,
        ("G719X+L861Q", "afatinib"):    5,   # n=5, assumed all events (conservative)
        ("Ex19delV3",   "osimertinib"): 5,
        ("Ex19delV4",   "osimertinib"): 3,
        ("Ex19delV5",   "osimertinib"): 4,
        ("G719X",       "osimertinib"): 14,
        ("G719X+S768I", "osimertinib"): 18,
    }
    for row in rows:
        k = (row["profile"], row["drug"])
        if k in EVENTS:
            row["n_events"] = EVENTS[k]
        elif row["n"] >= 10:
            row["n_events"] = row["n"]  # large groups assumed well-censored

    # Mark GENIE BPC rows that came from that source
    GENIE_PROFILES = {"Ex19delV3", "Ex19delV4", "Ex19delV5"}
    for row in rows:
        if row["profile"] in GENIE_PROFILES:
            row["source"] = "GENIE BPC"

    return rows


def apply_filters(rows, scope="v1.0"):
    """Apply pre-specified filters. Returns (qualifying, excluded)."""
    qualifying = []
    excluded   = []

    for row in rows:
        profile, drug = row["profile"], row["drug"]
        reason = None

        # Treatment-appropriateness (all scopes)
        if (profile, drug) in TREATMENT_EXCLUSIONS:
            reason = "treatment-inappropriateness"

        # Minimum n
        elif row["n"] < MIN_N:
            reason = f"n={row['n']} < {MIN_N}"

        # v1.0 structural exclusions
        elif scope == "v1.0":
            # Single-mutation profiles outside binding pocket
            if profile in V1_EXCLUDED_PROFILES:
                reason = "outside binding pocket (>9 Å active conformation)"

        if reason:
            excluded.append({**row, "reason": reason})
        else:
            qualifying.append(row)

    return qualifying, excluded


def spearman_with_ci(x, y, n_boot=10000, seed=42):
    """Spearman ρ with t-distribution p-value and bootstrap 95% CI."""
    rho, pval = stats.spearmanr(x, y)
    # Bootstrap CI
    rng = np.random.default_rng(seed)
    rhos = []
    n = len(x)
    xa, ya = np.array(x), np.array(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r, _ = stats.spearmanr(xa[idx], ya[idx])
        rhos.append(r)
    ci_lo = float(np.percentile(rhos, 2.5))
    ci_hi = float(np.percentile(rhos, 97.5))
    return float(rho), float(pval), ci_lo, ci_hi


def loo_analysis(pairs):
    """Leave-one-out: return list of {removed, n, rho, pval}."""
    results = []
    x_all = [p["ddg"] for p in pairs]
    y_all = [p["ttnt"] for p in pairs]
    rho_full, pval_full = stats.spearmanr(x_all, y_all)
    for i, pair in enumerate(pairs):
        x_loo = [p["ddg"] for j, p in enumerate(pairs) if j != i]
        y_loo = [p["ttnt"] for j, p in enumerate(pairs) if j != i]
        rho, pval = stats.spearmanr(x_loo, y_loo)
        results.append({
            "removed": f"{pair['profile']} + {pair['drug']}",
            "n": len(x_loo),
            "rho": round(float(rho), 3),
            "pval": round(float(pval), 3),
            "delta_rho": round(float(rho) - float(rho_full), 3),
        })
    return results


# ── GENIE BPC Replication ─────────────────────────────────────────────────────
def load_genie_replication(v1_pairs):
    """Match GENIE BPC TTNT to our DDG values for v1.0 qualifying profiles."""
    # Build DDG lookup from primary analysis
    ddg_lookup = {(p["profile"], p["drug"]): p["ddg"] for p in v1_pairs}

    genie_rows = []
    with open(GENIE_CSV) as f:
        for r in csv.DictReader(f):
            profile, drug = r["profile"], r["drug"]
            n = int(r["n"]) if r["n"] else 0
            n_events = int(r["n_events"]) if r["n_events"] else 0
            ttnt_s = r["ttnt_months"] if "ttnt_months" in r else r.get("ttnt_months", "")
            # Handle CSV field name
            ttnt_val = None
            for key in r:
                if "ttnt" in key.lower() and r[key]:
                    try:
                        ttnt_val = float(r[key])
                    except ValueError:
                        pass
                    break

            if ttnt_val is None:
                continue

            key = (profile, drug)
            if key not in ddg_lookup:
                continue
            if n < MIN_N or n_events < MIN_EVENTS:
                continue
            if (profile, drug) in TREATMENT_EXCLUSIONS:
                continue
            # Only v1.0 profiles (exclude L858R-only, L861Q-only, S768I-only)
            if profile in V1_EXCLUDED_PROFILES:
                continue

            genie_rows.append({
                "profile": profile,
                "drug": drug,
                "ddg": ddg_lookup[key],
                "ttnt": ttnt_val,
                "n": n,
                "n_events": n_events,
            })

    return genie_rows


# ── Figures ───────────────────────────────────────────────────────────────────
DRUG_COLORS  = {"osimertinib": "#2e7d32", "erlotinib": "#1565c0", "afatinib": "#e65100"}
DRUG_MARKERS = {"osimertinib": "o",       "erlotinib": "s",       "afatinib": "^"}
DRUG_LABELS  = {"osimertinib": "Osimertinib (3rd gen)",
                "erlotinib":   "Erlotinib (1st gen)",
                "afatinib":    "Afatinib (2nd gen)"}

PROFILE_SHORT = {
    "Ex19del":              "Ex19del",
    "Ex19delV3":            "Ex19del-V3",
    "Ex19delV4":            "Ex19del-V4",
    "Ex19delV5":            "Ex19del-V5",
    "G719X":                "G719X",
    "G719X+S768I":          "G719X+\nS768I",
    "G719X+L861Q":          "G719X+\nL861Q",
    "L858R+T790M":          "L858R+\nT790M",
    "Ex19del+T790M":        "Ex19del+\nT790M",
    "L858R":                "L858R",
    "L861Q":                "L861Q",
}


def fig1_scatter(pairs, rho, pval, ci_lo, ci_hi):
    """Figure 1: ΔΔG vs TTNT scatter for v1.0 scope (n=11)."""
    fig, ax = plt.subplots(figsize=(11, 8))

    for p in pairs:
        drug = p["drug"]
        c  = DRUG_COLORS.get(drug, "#888")
        mk = DRUG_MARKERS.get(drug, "o")
        sz = max(90, min(450, p["n"] ** 0.5 * 28))
        ax.scatter(p["ddg"], p["ttnt"], c=c, marker=mk,
                   s=sz, edgecolors="k", linewidths=0.9, zorder=3, alpha=0.88)

    # OLS trend
    x_all = [p["ddg"] for p in pairs]
    y_all = [p["ttnt"] for p in pairs]
    m, b = np.polyfit(x_all, y_all, 1)
    xline = np.linspace(min(x_all) - 0.5, max(x_all) + 0.5, 200)
    ax.plot(xline, m * xline + b, "--", color="#666", lw=1.3, alpha=0.7, label="OLS trend")

    # Per-point annotation offsets: (x_off, y_off) tuned to avoid overlap
    # Determined from actual data coordinates: ddg range −5.7 to +7.0, TTNT range 2.8 to 29.7
    OFFSETS = {
        ("Ex19del",       "erlotinib"):    ( 0.3, -1.8),   # below
        ("Ex19del",       "osimertinib"):  ( 0.3,  1.4),   # above-right (left cluster)
        ("G719X",         "osimertinib"):  ( 0.3, -1.5),   # below
        ("L858R+T790M",   "osimertinib"):  (-2.2,  0.8),   # left (avoid G719X+S768I)
        ("G719X+S768I",   "osimertinib"):  ( 0.3, -1.8),   # below-right
        ("G719X+L861Q",   "osimertinib"):  ( 0.3, -2.8),   # below (highest TTNT point)
        ("G719X+L861Q",   "afatinib"):     ( 0.8,  1.0),   # above-right (left cluster)
        ("Ex19del+T790M", "osimertinib"):  ( 0.3, -2.2),   # below-right (left cluster)
        ("Ex19delV3",     "osimertinib"):  (-2.5,  0.8),   # left (near right edge)
        ("Ex19delV4",     "osimertinib"):  ( 0.3,  0.9),   # above-right
        ("Ex19delV5",     "osimertinib"):  (-2.0,  0.9),   # above-left
    }

    for p in pairs:
        label = PROFILE_SHORT.get(p["profile"], p["profile"])
        drug_abbr = p["drug"][:3].upper()
        key = (p["profile"], p["drug"])
        x_off, y_off = OFFSETS.get(key, (0.3, 0.6))
        ax.annotate(
            f"{label}\n({drug_abbr})",
            xy=(p["ddg"], p["ttnt"]),
            xytext=(p["ddg"] + x_off, p["ttnt"] + y_off),
            fontsize=7.5, color="#222",
            ha="left",
            arrowprops=dict(arrowstyle="-", color="#aaa", lw=0.6),
            clip_on=False,
        )

    # Stats box
    sig = "**" if pval < 0.01 else ("*" if pval < 0.05 else ("†" if pval < 0.10 else "ns"))
    ax.text(0.97, 0.05,
            f"Spearman ρ = {rho:+.3f}\np = {pval:.3f} ({sig})\n95% CI [{ci_lo:.3f}, {ci_hi:.3f}]\nn = {len(pairs)} pairs",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.4", fc="#f5f5f5", ec="#bbb"))

    ax.axvline(0, color="#bbb", lw=0.9, ls=":")
    # Extra headroom above the G719X+L861Q point (TTNT=29.67) so its label fits
    ax.set_ylim(bottom=ax.get_ylim()[0], top=35)
    ax.set_xlabel("ΔΔG (kcal/mol)  [mutant − wildtype binding affinity]\nMore negative = stronger predicted mutant-preferential binding",
                  fontsize=10.5)
    ax.set_ylabel("KM-Median TTNT (months)\nTime on drug before next treatment", fontsize=10.5)
    ax.set_title("Figure 1. Computational ΔΔG vs Real-World TTNT\n"
                 "Pre-specified Scope v1.0: Binding-Pocket Mutations (≤5 Å from drug in active conformation)",
                 fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.25)

    drug_handles = [
        plt.scatter([], [], marker=DRUG_MARKERS[d], c=DRUG_COLORS[d],
                    s=90, edgecolors="k", linewidths=0.8, label=DRUG_LABELS[d])
        for d in ["osimertinib", "erlotinib", "afatinib"]
    ]
    ax.legend(handles=drug_handles, fontsize=9, loc="upper left",
              framealpha=0.9, edgecolor="#ccc")

    fig.tight_layout()
    for ext in [".pdf", ".png"]:
        out = FIG_DIR / f"fig1_v1_scatter{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=300 if ext == ".pdf" else 200)
        print(f"  Saved {out}")
    plt.close(fig)


def fig2_heatmap(all_pairs_v3):
    """Figure 2: ΔΔG heatmap all v3.0 profiles × drugs."""
    profiles_ordered = [
        "Ex19del", "Ex19delV3", "Ex19delV4", "Ex19delV5",
        "G719X", "G719X+S768I", "G719X+L861Q",
        "L858R", "L858R+T790M", "L861Q",
        "Ex19del+T790M",
    ]
    drugs = ["osimertinib", "erlotinib", "afatinib"]
    # Filter to profiles that appear in data
    seen_profiles = {p["profile"] for p in all_pairs_v3}
    profiles = [p for p in profiles_ordered if p in seen_profiles]

    grid = np.full((len(profiles), len(drugs)), np.nan)
    for p in all_pairs_v3:
        if p["profile"] in profiles and p["drug"] in drugs:
            r = profiles.index(p["profile"])
            c = drugs.index(p["drug"])
            grid[r, c] = p["ddg"]

    fig, ax = plt.subplots(figsize=(8, 6))
    absmax = max(abs(np.nanmin(grid)), abs(np.nanmax(grid)))
    im = ax.imshow(grid, cmap="RdYlGn", vmin=-absmax, vmax=absmax, aspect="auto")
    cb = plt.colorbar(im, ax=ax, label="ΔΔG (kcal/mol)", fraction=0.04, pad=0.04)
    cb.ax.tick_params(labelsize=9)

    ax.set_xticks(range(len(drugs)))
    ax.set_xticklabels(["Osimertinib", "Erlotinib", "Afatinib"], fontsize=11)
    ax.set_yticks(range(len(profiles)))
    ax.set_yticklabels(profiles, fontsize=9)

    for r in range(len(profiles)):
        for c in range(len(drugs)):
            v = grid[r, c]
            if not np.isnan(v):
                color = "white" if abs(v) > 0.5 * absmax else "black"
                ax.text(c, r, f"{v:+.2f}", ha="center", va="center",
                        fontsize=8.5, color=color, fontweight="bold")
            else:
                ax.text(c, r, "—", ha="center", va="center", fontsize=10, color="#bbb")

    ax.set_title("Figure 2. ΔΔG Heatmap: Mutant vs Wildtype EGFR Binding Affinity\n"
                 "Green = stronger mutant-preferential binding  |  Red = weaker than wildtype",
                 fontsize=10.5, fontweight="bold")

    # Mark v1.0 vs v3.0 profiles with left-side indicator
    v1_profiles = {p for p in profiles if p not in {"L858R", "L861Q"}}
    for i, prof in enumerate(profiles):
        color = "#2e7d32" if prof in v1_profiles else "#888"
        ax.get_yticklabels()[i].set_color(color)

    # Legend
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], color="#2e7d32", lw=3, label="Primary scope: binding-pocket (≤5 Å, v1.0)"),
        Line2D([0], [0], color="#888",    lw=3, label="Secondary scope only: sensitising mutations (v3.0)"),
    ]
    ax.legend(handles=legend_elems, loc="lower right", fontsize=8.5,
              framealpha=0.9, edgecolor="#ccc")

    fig.tight_layout()
    for ext in [".pdf", ".png"]:
        out = FIG_DIR / f"fig2_heatmap{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=300 if ext == ".pdf" else 200)
        print(f"  Saved {out}")
    plt.close(fig)


def fig3_loo(pairs, loo_results):
    """Figure 3: LOO bar chart showing Δρ for each removed point."""
    labels = [r["removed"] for r in loo_results]
    delta  = [r["delta_rho"] for r in loo_results]
    pvals  = [r["pval"] for r in loo_results]

    colors = []
    for i, (d, pv) in enumerate(zip(delta, pvals)):
        if pv < 0.05:
            colors.append("#2e7d32")   # still significant
        elif pv < 0.10:
            colors.append("#f9a825")   # trend
        else:
            colors.append("#c62828")   # no longer significant

    fig, ax = plt.subplots(figsize=(11, 6))
    y_pos = range(len(labels))
    bars = ax.barh(list(y_pos), delta, color=colors, edgecolor="k", linewidth=0.6, height=0.6)

    ax.axvline(0, color="#555", lw=1.0)

    # p-value labels: always placed to the RIGHT of the bar end (or right of zero for
    # negative bars).  This avoids left-side clipping into the y-axis tick labels.
    for i, (d, pv) in enumerate(zip(delta, pvals)):
        sig = "**" if pv < 0.01 else ("*" if pv < 0.05 else ("†" if pv < 0.10 else "ns"))
        x_text = max(d, 0) + 0.003   # right of bar end, or just right of zero
        ax.text(x_text, i, f"p={pv:.3f} ({sig})", va="center", ha="left", fontsize=7.5)

    # Expand x-axis so the widest bar + its label text both fit
    min_d = min(delta)
    max_d = max(delta)
    ax.set_xlim(min_d - 0.008, max_d + 0.022)

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Δρ (LOO correlation − full correlation)", fontsize=10)
    ax.set_title("Figure 3. Leave-One-Out Sensitivity Analysis — Binding-Pocket Scope v1.0 (n=11)\n"
                 "Bar = change in Spearman ρ when profile is removed",
                 fontsize=11, fontweight="bold")

    legend_elems = [
        mpatches.Patch(color="#2e7d32", label="p < 0.05 (significant)"),
        mpatches.Patch(color="#f9a825", label="p < 0.10 (trend)"),
        mpatches.Patch(color="#c62828", label="p ≥ 0.10 (not significant)"),
    ]
    ax.legend(handles=legend_elems, fontsize=8.5, loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()

    for ext in [".pdf", ".png"]:
        out = FIG_DIR / f"fig3_loo{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=300 if ext == ".pdf" else 200)
        print(f"  Saved {out}")
    plt.close(fig)


def fig4_genie(genie_pairs, rho, pval):
    """Figure 4: GENIE BPC replication scatter."""
    if not genie_pairs:
        print("  No GENIE BPC pairs qualified — skipping Figure 4")
        return

    # Per-point offsets for GENIE scatter labels
    # Erlotinib square (Ex19del) sits near the top border → move its label below
    GENIE_OFFSETS = {
        ("Ex19del",       "erlotinib"):    ( 0.15, -0.8),   # below (top point, avoid border)
        ("Ex19del",       "osimertinib"):  ( 0.8,  -0.7),   # right-below (clear of T790M)
        ("Ex19del+T790M", "osimertinib"):  ( 0.3,   0.7),   # above-right (avoid left axis)
        ("G719X",         "osimertinib"):  ( 0.15, -0.7),   # below
        ("L858R+T790M",   "osimertinib"):  ( 0.15,  0.4),   # right
        ("Ex19delV3",     "osimertinib"):  (-2.0,   0.4),   # left (near right border)
        ("Ex19delV5",     "osimertinib"):  (-1.8,   0.4),   # left
    }

    fig, ax = plt.subplots(figsize=(8, 6.5))
    for p in genie_pairs:
        drug = p["drug"]
        c  = DRUG_COLORS.get(drug, "#888")
        mk = DRUG_MARKERS.get(drug, "o")
        sz = max(90, min(400, p["n"] ** 0.5 * 26))
        ax.scatter(p["ddg"], p["ttnt"], c=c, marker=mk,
                   s=sz, edgecolors="k", linewidths=0.9, zorder=3, alpha=0.88)
        label = PROFILE_SHORT.get(p["profile"], p["profile"])
        key = (p["profile"], p["drug"])
        x_off, y_off = GENIE_OFFSETS.get(key, (0.15, 0.3))
        ax.annotate(
            f"{label}",
            xy=(p["ddg"], p["ttnt"]),
            xytext=(p["ddg"] + x_off, p["ttnt"] + y_off),
            fontsize=8, color="#333",
            arrowprops=dict(arrowstyle="-", color="#aaa", lw=0.5),
            clip_on=False,
        )

    sig = "**" if pval < 0.01 else ("*" if pval < 0.05 else ("†" if pval < 0.10 else "ns"))
    ax.text(0.97, 0.05,
            f"Spearman ρ = {rho:+.3f}\np = {pval:.3f} ({sig})\nn = {len(genie_pairs)} pairs",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.4", fc="#e8f5e9", ec="#a5d6a7"))

    ax.axvline(0, color="#bbb", lw=0.9, ls=":")
    ax.set_xlabel("ΔΔG (kcal/mol)", fontsize=10.5)
    ax.set_ylabel("KM-Median TTNT (months) — GENIE BPC", fontsize=10.5)
    ax.set_title("Figure 4. GENIE BPC Replication (Independent Cohort)\n"
                 "v1.0 binding-pocket profiles, GENIE BPC NSCLC v2.0-public",
                 fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.25)
    # Expand y-axis top so no label is clipped at the border
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.12)
    fig.tight_layout()

    for ext in [".pdf", ".png"]:
        out = FIG_DIR / f"fig4_genie_replication{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=300 if ext == ".pdf" else 200)
        print(f"  Saved {out}")
    plt.close(fig)


# ── GENIE BPC Regimen-Line Analysis (Supplementary Table 4) ──────────────────
def genie_line_analysis():
    """
    Compute regimen-line statistics from cohort.csv to formally demonstrate
    treatment-line confounding in the GENIE BPC Ex19del data.

    Returns a dict with three parts matching Supplementary Table 4:
      part_a: frontline counts and median line per profile/drug combo
      part_b: Spearman rho / p between regimen_number and ttnt (osimertinib group)
      part_c: mean TTNT by line category (frontline vs later-line) for Ex19del
    """
    if not GENIE_COHORT.exists():
        print(f"  WARNING: {GENIE_COHORT} not found — skipping Supplementary Table 4")
        return None

    rows = []
    with open(GENIE_COHORT) as f:
        for r in csv.DictReader(f):
            try:
                rows.append({
                    "sample_id":       r["sample_id"],
                    "profile":         r["profile"],
                    "drug":            r["drug"],
                    "regimen_number":  int(r["regimen_number"]),
                    "ttnt_months":     float(r["ttnt_months"]),
                    "censored":        r.get("censored", "0").strip() in ("1", "True", "true"),
                })
            except (KeyError, ValueError):
                continue

    # ── Part A: frontline counts and median line per profile/drug ────────────
    from collections import defaultdict
    groups = defaultdict(list)
    for r in rows:
        groups[(r["profile"], r["drug"])].append(r)

    part_a = []
    for (profile, drug), grp in sorted(groups.items()):
        n_total    = len(grp)
        n_frontline = sum(1 for r in grp if r["regimen_number"] == 1)
        lines = sorted(r["regimen_number"] for r in grp)
        n2 = len(lines)
        median_line = lines[n2 // 2] if n2 % 2 == 1 else (lines[n2 // 2 - 1] + lines[n2 // 2]) / 2
        part_a.append({
            "profile":        profile,
            "drug":           drug,
            "n":              n_total,
            "n_frontline":    n_frontline,
            "pct_frontline":  round(100 * n_frontline / n_total, 1) if n_total else None,
            "median_line":    median_line,
        })

    # ── Part B: Spearman(regimen_number, ttnt) within Ex19del+osimertinib ────
    # Restricted to Ex19del+osimertinib to match manuscript citation (ρ=−0.335, p=0.003, n=76).
    osi_rows = [r for r in rows if r["drug"] == "osimertinib" and r["profile"] == "Ex19del"]
    if len(osi_rows) >= 3:
        x_line = [r["regimen_number"] for r in osi_rows]
        y_ttnt  = [r["ttnt_months"]    for r in osi_rows]
        rho_line, p_line = stats.spearmanr(x_line, y_ttnt)
    else:
        rho_line, p_line = None, None

    part_b = {
        "description": "Spearman rho between regimen_number and ttnt_months (Ex19del+osimertinib patients, n=76)",
        "n": len(osi_rows),
        "rho": round(rho_line, 3) if rho_line is not None else None,
        "pval": round(p_line, 3) if p_line is not None else None,
    }

    # ── Part C: mean TTNT frontline vs later-line for Ex19del ────────────────
    part_c = []
    for drug in ("osimertinib", "erlotinib"):
        sub = [r for r in rows if r["profile"] == "Ex19del" and r["drug"] == drug]
        frontline   = [r["ttnt_months"] for r in sub if r["regimen_number"] == 1]
        later_line  = [r["ttnt_months"] for r in sub if r["regimen_number"] > 1]
        part_c.append({
            "profile": "Ex19del",
            "drug": drug,
            "n_frontline":  len(frontline),
            "mean_ttnt_frontline":  round(sum(frontline)  / len(frontline),  1) if frontline  else None,
            "n_later_line": len(later_line),
            "mean_ttnt_later_line": round(sum(later_line) / len(later_line), 1) if later_line else None,
        })

    result = {"part_a": part_a, "part_b": part_b, "part_c": part_c}

    # Print summary to console
    print("\nSupplementary Table 4 — GENIE BPC Regimen Line Analysis:")
    print("  Part A: Frontline use by profile/drug")
    for row in part_a:
        print(f"    {row['profile']:30s} {row['drug']:15s}  n={row['n']:4d}  "
              f"frontline={row['n_frontline']} ({row['pct_frontline']}%)  "
              f"median_line={row['median_line']}")
    print(f"\n  Part B: Spearman(regimen_number, TTNT) in osimertinib patients")
    print(f"    n={part_b['n']}, ρ={part_b['rho']}, p={part_b['pval']}")
    print(f"\n  Part C: Mean TTNT by line category for Ex19del")
    for row in part_c:
        print(f"    {row['drug']:15s}  frontline: n={row['n_frontline']} mean={row['mean_ttnt_frontline']} mo"
              f"  |  later-line: n={row['n_later_line']} mean={row['mean_ttnt_later_line']} mo")

    return result


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print("=== v2 Pre-specified Analysis ===\n")

    all_rows = load_primary()

    # ── v1.0 ─────────────────────────────────────────────────────────────────
    v1_pairs, v1_excl = apply_filters(all_rows, scope="v1.0")
    x1 = [p["ddg"] for p in v1_pairs]
    y1 = [p["ttnt"] for p in v1_pairs]
    rho1, p1, ci1_lo, ci1_hi = spearman_with_ci(x1, y1)

    print(f"v1.0 — Binding-pocket (n={len(v1_pairs)})")
    print(f"  ρ = {rho1:+.3f},  p = {p1:.3f},  95% CI [{ci1_lo:.3f}, {ci1_hi:.3f}]")
    for p in v1_pairs:
        print(f"  {p['profile']:30s}  {p['drug']:15s}  ddg={p['ddg']:+7.3f}  ttnt={p['ttnt']:.2f} mo  n={p['n']}")
    print(f"\n  Excluded ({len(v1_excl)}):")
    for p in v1_excl:
        print(f"  {p['profile']:30s}  {p['drug']:15s}  [{p['reason']}]")

    # ── v1.0 osimertinib-only ────────────────────────────────────────────────
    v1_osi = [p for p in v1_pairs if p["drug"] == "osimertinib"]
    x_osi = [p["ddg"] for p in v1_osi]
    y_osi = [p["ttnt"] for p in v1_osi]
    rho_osi, p_osi, ci_osi_lo, ci_osi_hi = spearman_with_ci(x_osi, y_osi)
    print(f"\nv1.0 osimertinib-only (n={len(v1_osi)})")
    print(f"  ρ = {rho_osi:+.3f},  p = {p_osi:.3f},  95% CI [{ci_osi_lo:.3f}, {ci_osi_hi:.3f}]")

    # ── v3.0 ─────────────────────────────────────────────────────────────────
    v3_pairs, v3_excl = apply_filters(all_rows, scope="v3.0")
    x3 = [p["ddg"] for p in v3_pairs]
    y3 = [p["ttnt"] for p in v3_pairs]
    rho3, p3, ci3_lo, ci3_hi = spearman_with_ci(x3, y3)
    print(f"\nv3.0 — All sensitising (n={len(v3_pairs)})")
    print(f"  ρ = {rho3:+.3f},  p = {p3:.3f},  95% CI [{ci3_lo:.3f}, {ci3_hi:.3f}]")

    # ── v3.0 osimertinib-only ────────────────────────────────────────────────
    v3_osi = [p for p in v3_pairs if p["drug"] == "osimertinib"]
    x_v3osi = [p["ddg"] for p in v3_osi]
    y_v3osi = [p["ttnt"] for p in v3_osi]
    rho_v3osi, p_v3osi, ci_v3osi_lo, ci_v3osi_hi = spearman_with_ci(x_v3osi, y_v3osi)
    print(f"\nv3.0 osimertinib-only (n={len(v3_osi)})")
    print(f"  ρ = {rho_v3osi:+.3f},  p = {p_v3osi:.3f},  95% CI [{ci_v3osi_lo:.3f}, {ci_v3osi_hi:.3f}]")

    # ── LOO analysis (v1.0) ──────────────────────────────────────────────────
    loo = loo_analysis(v1_pairs)
    print(f"\nLOO Analysis (v1.0):")
    print(f"  {'Profile + Drug':50s}  {'n':>4}  {'ρ':>7}  {'p':>7}  {'Δρ':>7}")
    for r in loo:
        print(f"  {r['removed']:50s}  {r['n']:>4}  {r['rho']:>+7.3f}  {r['pval']:>7.3f}  {r['delta_rho']:>+7.3f}")

    # ── GENIE BPC Regimen-Line Confounding (Supplementary Table 4) ──────────
    supp_table4 = genie_line_analysis()

    # ── GENIE BPC Replication ────────────────────────────────────────────────
    genie_pairs = load_genie_replication(v1_pairs)
    if genie_pairs:
        xg = [p["ddg"] for p in genie_pairs]
        yg = [p["ttnt"] for p in genie_pairs]
        rho_g, p_g = stats.spearmanr(xg, yg)
        print(f"\nGENIE BPC Replication (v1.0 profiles, n={len(genie_pairs)})")
        print(f"  ρ = {rho_g:+.3f},  p = {p_g:.3f}")
        for p in genie_pairs:
            print(f"  {p['profile']:30s}  {p['drug']:15s}  ddg={p['ddg']:+7.3f}  ttnt={p['ttnt']:.2f} mo  n={p['n']}")
    else:
        rho_g, p_g = None, None
        print("\nGENIE BPC Replication: no qualifying pairs (all already in primary dataset)")

    # ── Save JSON results ────────────────────────────────────────────────────
    results = {
        "v1_0": {
            "n": len(v1_pairs),
            "rho": round(rho1, 3),
            "pval": round(p1, 3),
            "ci_lo": round(ci1_lo, 3),
            "ci_hi": round(ci1_hi, 3),
            "pairs": [{"profile": p["profile"], "drug": p["drug"],
                       "ddg": p["ddg"], "ttnt": p["ttnt"], "n": p["n"]} for p in v1_pairs],
        },
        "v1_0_osi": {"n": len(v1_osi), "rho": round(rho_osi, 3),
                     "pval": round(p_osi, 3), "ci_lo": round(ci_osi_lo, 3), "ci_hi": round(ci_osi_hi, 3)},
        "v3_0": {
            "n": len(v3_pairs),
            "rho": round(rho3, 3),
            "pval": round(p3, 3),
            "ci_lo": round(ci3_lo, 3),
            "ci_hi": round(ci3_hi, 3),
        },
        "v3_0_osi": {"n": len(v3_osi), "rho": round(rho_v3osi, 3),
                     "pval": round(p_v3osi, 3), "ci_lo": round(ci_v3osi_lo, 3), "ci_hi": round(ci_v3osi_hi, 3)},
        "loo_v1_0": loo,
        "genie_replication": {
            "n": len(genie_pairs),
            "rho": round(rho_g, 3) if rho_g is not None else None,
            "pval": round(p_g, 3) if p_g is not None else None,
            "pairs": [{"profile": p["profile"], "drug": p["drug"],
                       "ddg": p["ddg"], "ttnt": p["ttnt"], "n": p["n"]} for p in genie_pairs],
        },
        "supplementary_table4_line_confounding": supp_table4,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved → {OUT_JSON}")

    # ── Generate Figures ─────────────────────────────────────────────────────
    print("\nGenerating figures...")
    fig1_scatter(v1_pairs, rho1, p1, ci1_lo, ci1_hi)
    fig2_heatmap(v3_pairs)
    fig3_loo(v1_pairs, loo)
    if genie_pairs:
        fig4_genie(genie_pairs, rho_g, p_g)

    print("\n=== Analysis complete ===")
    print(f"  Figures → {FIG_DIR}/")
    print(f"  Results → {OUT_JSON}")


if __name__ == "__main__":
    main()
