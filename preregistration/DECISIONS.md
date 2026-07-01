# Study 2 — Analytical Decisions Log

**Date:** 2026-06-26  
**Author:** Ahmad Abboud, Algentics Research  
**Purpose:** Single reference document for every analytical decision in this study —  
what was chosen, what was rejected, and why. Written before v2 analysis runs.  
This document is the source of truth for the Methods section when writing the paper.

Anyone reading this should be able to understand every choice without needing  
to read any other file.

---

## Decision 1 — Clinical outcome: TTNT, not PFS or OS

**Chosen:** Time to Next Treatment (TTNT)  
**Definition:** Days from start of a given treatment line to start of the next treatment line (or death). Converted to months at 30.4375 days/month. Patients without a next line are censored at last recorded clinical event.

**Why not PFS (progression-free survival):**  
PFS requires a radiological progression date — a CT scan read as "worsening" with a documented date. The public MSK-CHORD cBioPortal API exposes `HAS_CANCER` (cancer present: yes/no), not radiological progression. In active cancer patients, every scan returns `HAS_CANCER = Y`. When PFS was tested using this field, 98% of observations appeared as early events, producing median "PFS" of 0.6–1.8 months across all profiles — a data artifact, not a real signal. PFS is therefore not computable from this data source.

**Why not OS (overall survival):**  
OS measures time until death. It is driven by many factors beyond whether the first TKI worked — second and third line therapies, performance status, comorbidities. It would dilute the specific signal we are measuring: "how well did this drug bind and suppress this mutation."

**Why TTNT is valid here:**  
1. For NSCLC patients on TKIs, the interval between progression and next treatment start is short — clinicians switch quickly. TTNT ≈ PFS in this setting.  
2. Empirical validation: osimertinib TTNT is 2.2× longer than erlotinib for the same profiles (Ex19del: 11.20 vs 4.83 mo; L858R: 10.68 vs 4.89 mo), consistent with the FLAURA trial direction. TTNT is capturing real treatment efficacy differences.  
3. TTNT is a validated PFS surrogate in real-world NSCLC registry data (Griffith et al. 2019, Adv Ther).

**Known limitation:**  
A patient who stopped treatment due to toxicity — not because the drug stopped working — also contributes a short TTNT. This adds noise. It is disclosed as Limitation 2 in the manuscript.

---

## Decision 2 — Minimum patient threshold: n ≥ 5 per pair

**Chosen:** n ≥ 5 patients per profile × drug pair with a computable KM median.

**Why 5 and not lower:**  
A KM median requires S(t) to cross 0.50. With n = 2–3 patients and heavy censoring, the KM curve may never cross 0.50, making the median undefined. Even when computable, a median from 2–3 patients is statistically meaningless and clinically unreliable — one outlier patient changes the median dramatically. n = 5 is the minimum at which a KM median has any interpretive value.

**Why 5 and not higher (e.g., 10):**  
This is a proof-of-concept study. Rare profiles (G719X compounds, C797S triples) are inherently small. Setting n ≥ 10 would exclude exactly the profiles that are most scientifically interesting — the compound mutations with no clinical trial data. n = 5 is a deliberate trade-off between reliability and coverage.

**What this excludes:**  
- G719X + Afatinib (n=4 in MSK-CHORD) — below threshold  
- G719X+S768I + Afatinib (n=3) — below threshold  
- C797S+L858R+T790M + Osimertinib (n=2) — below threshold  
- T790M alone (n=1) — below threshold and biologically incorrect (T790M is always acquired, never de novo driver)

---

## Decision 3 — Drug selection: osimertinib, erlotinib, afatinib only

**Chosen:** Osimertinib (3rd gen), Erlotinib (1st gen), Afatinib (2nd gen).

**Why these three:**  
These are the three EGFR TKIs with sufficient patient numbers in MSK-CHORD and GENIE BPC, and for which validated co-crystal structures exist for wildtype EGFR (needed for ΔΔG computation):
- Osimertinib: PDB 4ZAU (WT EGFR + AZD9291)
- Erlotinib: PDB 4HJO (WT EGFR + erlotinib analog)
- Afatinib: PDB 4G5J (WT EGFR + BIBW2992)

**Why not gefitinib:**  
Gefitinib is common in East Asian cohorts but rare in MSK-CHORD (US/Western registry). Insufficient patient numbers for any profile to meet n ≥ 5. No co-crystal structure used as WT reference in this pipeline.

**Why not dacomitinib or other 2nd/3rd gen TKIs:**  
Same reason — insufficient patient counts in both databases.

**Important constraint — ΔΔG not comparable across drugs:**  
Each drug uses a different WT reference structure. ΔG(wildtype) differs per drug (−4.10 kcal/mol for osimertinib, −6.95 for erlotinib, −1.47 for afatinib). ΔΔG values can only be compared within a single drug. Mixing drugs in one correlation adds noise and must be interpreted carefully.

---

## Decision 4 — Treatment-appropriateness filter

**Chosen:** Exclude profile × drug pairs that are clinically inappropriate.

**What is excluded and why:**

| Excluded pair | Reason |
|---|---|
| T790M-containing profile + Erlotinib | Erlotinib IC50 is ~100× worse against T790M. Clinically ineffective — any TTNT reflects toxicity-driven discontinuation, not drug efficacy. |
| T790M-containing profile + Afatinib | Afatinib does not adequately cover T790M at clinical doses. Same problem. |
| Ex19del+C797S + Osimertinib | C797S is the osimertinib resistance mutation. Treating this profile with osimertinib is clinically irrational. |

**Why this is not data-driven exclusion:**  
These exclusions are based on known pharmacology (IC50 data, resistance mechanisms) established before any TTNT data was seen. They would be listed in any study protocol before data collection.

---

## Decision 5 — Wildtype reference structures for ΔΔG

**Chosen:** Drug-matched co-crystal structures.

| Drug | PDB ID | Co-crystal ligand | Resolution |
|---|---|---|---|
| Osimertinib | 4ZAU | YY3 (AZD9291) | 2.5 Å |
| Erlotinib | 4HJO | AQ4 (erlotinib analog) | 2.7 Å |
| Afatinib | 4G5J | 0WN (BIBW2992) | 2.8 Å |

**Why drug-matched:**  
Each TKI binds EGFR in a slightly different conformation. Using a drug-matched co-crystal structure ensures the wildtype receptor geometry is appropriate for the docking autobox of that specific drug. Using a single WT structure for all three drugs would bias the autobox toward one drug's binding mode.

**Special case — Ex19del+T790M:**  
The Boltz-2-predicted structure for this profile failed quality control (active-site coordinate misalignment). PDB 6LUD was used as a backbone template: the exon 19 deletion was applied, T790M was introduced, and C797S was reverted to CYS to enable covalent docking. 6LUD is NOT a wildtype reference — it is a structure-building template for this one profile only.

---

## Decision 6 — Statistical method: Spearman rank correlation

**Chosen:** Spearman ρ with t-distribution p-value (scipy.stats.spearmanr, df = n−2).

**Why Spearman and not Pearson:**  
1. The ΔΔG–TTNT relationship may be non-linear — Spearman captures monotonic association without assuming linearity.  
2. KM median TTNT is an ordinal summary — it ranks profiles but does not have the distributional properties required for Pearson.  
3. Outlier profiles (e.g., G719X+L861Q with TTNT = 29.67 months from 3 events) would disproportionately influence a parametric test.

**Why t-distribution and not z-approximation:**  
For small n (n < 30), the normal (z) approximation for Spearman p-values is anti-conservative — it gives p-values that are too small. The t-distribution with df = n−2 is correct for small samples and is the method used by scipy.stats.spearmanr. This was identified as a bug in v1 exploratory: the custom implementation used z, giving p = 0.030; the correct t-distribution gives p = 0.055 for the same ρ = −0.566 at n = 12.

**Bootstrap 95% CI:**  
10,000 resamples with replacement. 2.5th and 97.5th percentiles. Used because the theoretical CI for Spearman ρ assumes normality which does not hold for small n.

**Significance threshold:**  
Two-tailed α = 0.05. Pre-specified. Not adjusted post-hoc.

---

## Decision 7 — Unit of analysis: profile × drug pair, not individual patient

**Chosen:** Each (mutation profile, drug) combination is one data point. The y-value is the KM median TTNT across all patients in that group.

**Why not patient-level regression:**  
Patient-level analysis would require adjusting for clinical confounders (age, stage, ECOG, prior therapies, co-mutations). This information is either unavailable or unreliable in the public registry data. Profile-level KM medians smooth over individual patient variation and give a stable estimate of the "typical" outcome for that profile × drug combination.

**Known limitation:**  
Treating profile medians as independent data points is appropriate for correlation but prevents patient-level confounder adjustment. This is disclosed in the manuscript as Limitation 7.

---

## Decision 8 — Docking parameters

**Chosen:**
- GNINA v1.3.2
- Scoring: AutoDock Vina (`--cnn_scoring none`) — not CNN scoring
- Exhaustiveness: 32
- Autobox: 6 Å padding around reference ligand
- Replicates: 3 seeds per pair (42, 43, 44); report mean
- Covalent docking: for osimertinib and afatinib against C797-intact profiles

**Why AutoDock Vina scoring and not CNN:**  
CNN (convolutional neural network) scoring in GNINA is trained on PDBbind data and can overfit to training-set binding modes. For novel mutant structures generated by Boltz-2 (not in any training set), the physics-based Vina score is more reliable and reproducible.

**Why 3 seeds:**  
Docking is stochastic. 3 replicates with different random seeds quantify reproducibility. Mean SD across replicates was 0.14 kcal/mol in v1 — confirming the results are stable.

**Why covalent docking for osimertinib and afatinib:**  
Both drugs form a covalent bond with C797 via Michael addition of an acrylamide warhead. Non-covalent docking ignores this bond and systematically underestimates binding energy for these drugs. Profiles with C797S (the cysteine is mutated to serine — no thiol) receive non-covalent docking automatically.

---

## Decision 9 — Minimum events for KM median reliability

**Proposed for v2 (not applied in v1):**  
Exclude any profile × drug pair where n_events < 4.

**Rationale:**  
A KM median with fewer than 4 events is highly unstable — a single censored-vs-event reclassification changes the median substantially. G719X+L861Q osimertinib (3 events, n=6) is the clearest example: median TTNT = 29.67 months with CI not estimable. This data point is the most influential in the LOO analysis and likely reflects censoring inflation, not true treatment durability.

**Impact:**  
Under this rule, G719X+L861Q osimertinib would be excluded from the primary analysis. G719X+L861Q afatinib (n=5, events unknown from current data — needs verification) may also be excluded.

**Why this improves the analysis:**  
Removing unreliable data points reduces noise without cherry-picking — the rule applies mechanically to all pairs, not just the one that happens to be an outlier.

---

## Summary table — all decisions at a glance

| # | Decision | Chosen | Key reason |
|---|---|---|---|
| 1 | Clinical outcome | TTNT | PFS not computable from public API; TTNT is validated surrogate |
| 2 | Minimum n | ≥ 5 patients | KM median not reliable below this |
| 3 | Drugs | Osimertinib, erlotinib, afatinib | Only drugs with sufficient n and validated WT crystal structures |
| 4 | Treatment-appropriateness | Exclude clinically irrational pairs | Based on known pharmacology, not data |
| 5 | WT reference | Drug-matched co-crystal | Ensures correct binding mode geometry per drug |
| 6 | Statistics | Spearman + t-distribution + bootstrap CI | Correct for small n, non-linear, ordinal data |
| 7 | Unit of analysis | Profile × drug pair (KM median) | Patient-level confounders unavailable |
| 8 | Docking | GNINA v1.3.2, Vina scoring, 3 seeds | Physics-based, reproducible, covalent where applicable |
| 9 | Minimum events | ≥ 4 events per KM median (proposed) | Removes unstable KM estimates; applied in v2, not v1 |
