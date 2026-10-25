# Computational EGFR TKI Binding Affinity Correlates With Real-World TTNT

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Pre-registered](https://img.shields.io/badge/Analysis-Pre--specified-green)](preregistration/SCOPE.md)

**Abboud A, Dalloul Z, Dalloul I, Abdelsalam M.** *Computational Binding Affinity Correlates With Real-World Time to Next Treatment Across EGFR TKI–Mutation Profiles: A Pre-specified Correlation Study.* Algentics Research, 2026.

---

## What this study shows

We tested whether the computational tightness of a drug binding to a mutant EGFR protein
(ΔΔG, kcal/mol) predicts how long patients stayed on that drug before switching
(time to next treatment, TTNT, months) — across 11 mutation–drug profiles defined
by a pre-specified structural criterion.

**Primary result (v1.0 binding-pocket scope, n = 11 pairs):**

> ρ = −0.673, p = 0.023, bootstrap 95% CI: −0.886 to −0.174

More negative ΔΔG (stronger predicted binding) → longer real-world treatment duration.

**Secondary result (v3.0 all-sensitising scope, n = 15 pairs):**

> ρ = −0.643, p = 0.010, bootstrap 95% CI: −0.841 to −0.251

**Biochemical validation (ΔΔG vs ChEMBL pIC50, n = 5 profiles):**

> ρ = +0.900, p = 0.037 — docking scores align with enzymatic binding assays.

**GENIE BPC independent replication (n = 7 pairs, 5 fully independent):** ρ = +0.143, p = 0.760,
bootstrap 95% CI: −0.750 to +0.750 — non-confirming, but with a CI spanning 1.50 ρ units this
attempt is statistically uninformative at n = 7. The manuscript attributes the non-replication to
treatment-line enrollment bias in the GENIE BPC cohort (§3.5 and Supplementary Table 4).

**Patient-level Cox model (n = 576 treatment lines, 404 events):** ΔΔG hazard ratio 1.16
(95% CI 1.05–1.28, p = 0.0031) after adjustment for line of therapy, age, sex and stage —
a ~6% change from the unadjusted HR of 1.09, supporting the primary association at the level of
individual treatment lines rather than profile aggregates.

Leave-one-out analysis: all 11 subsets maintain the negative direction; 10/11 remain significant at p < 0.05.

---

## Repository structure

```
egfr-tki-binding-affinity-ttnt/
├── README.md                          ← this file
├── LICENSE                            ← MIT
├── requirements.txt                   ← Python dependencies
│
├── analysis/
│   ├── run_analysis.py                ← single script that reproduces all results + figures
│   └── provenance/                    ← record of two steps that need non-shipped inputs
│       ├── README.md                  ← what they did and what they require
│       ├── redock_erlotinib_1M17.py   ← 1M17 active-conformation erlotinib re-dock
│       └── run_patientlevel_cox.py    ← §3.7 patient-level Cox models
│
├── preregistration/                   ← written and dated BEFORE any analysis (2026-06-26)
│   ├── SCOPE.md                       ← structural inclusion criteria (v1.0–v3.0 decision tree)
│   ├── DECISIONS.md                   ← 9 analytical decisions with rationale
│   └── distance_measurements.md       ← measured atomic distances justifying v1.0 scope
│
├── data/
│   ├── README.md                      ← data sources and access instructions
│   ├── ddg_values.csv                 ← ΔΔG per profile × drug with scope annotation
│   ├── ddg_summary.csv                ← per-seed docking ΔG means and SDs
│   ├── merged_ddg_ttnt.csv            ← ΔΔG + MSK-CHORD TTNT, profile-level (shareable)
│   └── genie_bpc_profile_summary.csv  ← GENIE BPC TTNT, profile-level (shareable)
│
├── results/
│   ├── all_results.json               ← machine-readable results (all analyses)
│   ├── correlation_results.md         ← human-readable results summary
│   └── patientlevel_cox_results.json  ← §3.7 Cox model output
│
└── figures/
    ├── fig1_v1_scatter.png/.pdf       ← Figure 1: ΔΔG vs TTNT scatter (primary, n=11)
    ├── fig2_heatmap.png/.pdf          ← Figure 2: ΔΔG heatmap all profiles × drugs
    ├── fig3_loo.png/.pdf              ← Figure 3: leave-one-out sensitivity analysis
    └── fig4_genie_replication.png/.pdf ← Figure 4: GENIE BPC replication scatter
```

---

## How to reproduce

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

Requires Python ≥ 3.10. Dependencies: `numpy`, `matplotlib`, `scipy`.

### 2. Run the analysis

```bash
python analysis/run_analysis.py
```

This regenerates all four figures and `results/all_results.json` from the data files
in `data/`. Runtime: under 30 seconds on any modern laptop.

### 3. What this repository does not include

Two patient-level tables are deliberately not shipped, because both source datasets
are licensed CC BY-NC-ND and a derived patient-level table would be a redistribution:

- the **GENIE BPC cohort file** needed for Supplementary Table 4 (Synapse data access
  agreement; `analysis/run_analysis.py` detects its absence and skips that step);
- the **MSK-CHORD per-line table** needed to refit the §3.7 Cox models. The fitted
  model output is shipped in `results/patientlevel_cox_results.json`, and the table
  is reconstructible from cBioPortal study `msk_chord_2024` without registration.

Everything else reproduces from a clean checkout. See `analysis/provenance/README.md`.

---

## Computational pipeline

| Step | Tool | Key parameters |
|---|---|---|
| Structure prediction | Boltz-2 (AlphaFold3-class, Neurosnap API) | UniProt P00533, residues 696–1022 |
| Structure preparation | PDBFixer + OpenMM | AMBER14, 500 steps, implicit solvent, pH 7.4 |
| Molecular docking | GNINA v1.3.2 | `--cnn_scoring none --exhaustiveness 32`; covalent for osimertinib/afatinib on C797-intact profiles |
| Wildtype references | 4ZAU (osimertinib, active), 1M17 (erlotinib, active), 4G5J (afatinib, active) | Drug-matched |
| Replicates | 3 seeds (42, 43, 44) | Report mean ΔG |
| Statistics | Spearman rank | `scipy.stats.spearmanr`; bootstrap CI 10,000 resamples, seed 42 |

ΔΔG = ΔG(mutant) − ΔG(wildtype). More negative = stronger mutant-preferential binding.
ΔΔG values are not comparable across drugs (each uses a different wildtype reference structure).

---

## Pre-specification

The structural inclusion criterion (scope v1.0), all analytical decisions, and the
secondary scope (v3.0) were written and committed before any correlation was computed.
See `preregistration/` for the full documentation, dated 2026-06-26.

The analysis plan defined four scope versions as a sequential decision tree:
v1.0 (5 Å direct contact, primary) → v1.1 (8 Å fallback) → v2.0 (structural element fallback)
→ v3.0 (all OncoKB sensitising mutations, secondary). Because v1.0 reached significance
(p = 0.023), v1.1 and v2.0 were not triggered. Only v1.0 and v3.0 appear in the results.

---

## Data sources

Clinical data are from publicly available de-identified datasets:

- **MSK-CHORD 2024** — cBioPortal (study ID: `msk_chord_2024`), no registration required.
- **GENIE BPC NSCLC v2.0-public** — Synapse accession `syn27056697`, free registration required.

Neither dataset contains patient names or direct identifiers. See `data/README.md`.

---

## Citation

```
Abboud A, Dalloul Z, Dalloul I, Abdelsalam M. Computational Binding Affinity
Correlates With Real-World Time to Next Treatment Across EGFR TKI–Mutation
Profiles: A Pre-specified Correlation Study. Algentics Research, 2026.
```

---

## License

Code is released under the [MIT License](LICENSE).  
Pre-registration documents (`preregistration/`) are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
