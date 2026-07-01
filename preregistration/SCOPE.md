# Study 2 — Pre-specified Scope Definition

**Date locked:** 2026-06-26  
**Author:** Ahmad Abboud, Algentics Research  
**Rule:** This file is written and dated BEFORE any correlation analysis is run.  
Scope version must be incremented and re-dated if the inclusion criteria change.  
The version that produced the final reported result must be cited in the Methods section.

---

## Hypothesis (unchanged across all scope versions)

> Computational binding affinity (ΔΔG, kcal/mol) of an EGFR TKI to a patient-specific mutation  
> profile is negatively correlated with real-world time to next treatment (TTNT):  
> more favorable binding (more negative ΔΔG) → longer treatment duration.

---

## What does NOT change regardless of scope version

These rules apply to every version below:

| Rule | Detail |
|---|---|
| **ΔΔG definition** | ΔG(mutant) − ΔG(wildtype), kcal/mol. More negative = stronger mutant-preferential binding. |
| **TTNT definition** | Kaplan-Meier median time from treatment start to next treatment start (or death). Not a survival measure. |
| **Minimum n** | n ≥ 5 patients per profile × drug pair with a computable KM median. |
| **Data sources** | MSK-CHORD 2024 (cBioPortal msk_chord_2024) and GENIE BPC NSCLC v2.0-public (Synapse syn27056697). All qualifying pairs from both sources are included. No selective inclusion. |
| **Drugs** | Osimertinib, erlotinib, afatinib. No other drugs (insufficient docking validation). |
| **Statistics** | Spearman rank correlation. p-values via scipy.stats.spearmanr (t-distribution, df = n−2). Bootstrap 95% CI: 10,000 resamples. Two-tailed α = 0.05. |
| **One primary analysis** | One correlation per scope version. No post-hoc subsetting after seeing the result. |

---

## Scope versions

### v1.0 — Strict direct contact (5 Å criterion)

**Definition:**  
Include a mutation profile if ALL mutations in the profile satisfy:  
> The mutated residue's Cα atom is within **5 Å** of the bound ligand heavy atoms in the cognate drug-matched crystal structure (4ZAU for osimertinib, 4HJO for erlotinib, 4G5J for afatinib).

For compound profiles (e.g., Ex19del+T790M), at least one mutation must satisfy the criterion.  
Deletion mutations (Ex19del variants) qualify if any residue in the deleted segment falls within 5 Å.

**Structural basis:**  
5 Å is the standard cutoff for direct binding-site contact in structural biology.  
Residues beyond this distance do not make van der Waals or electrostatic contact with the drug.

**Expected inclusions under v1.0:**

| Profile | Key residue(s) | Within 5 Å of ligand? | Include? |
|---|---|---|---|
| Ex19del (all variants) | E746–A750 deletion (β3-αC loop) | Yes — loop lines the pocket entrance | ✓ |
| L858R | L858 (activation loop) | ~10–12 Å — does NOT directly contact ligand | ✗ |
| T790M | T790 (gatekeeper) | Yes — gatekeeper directly contacts drug | ✓ (via compound profiles) |
| G719X | G719 (P-loop) | Yes — P-loop contacts drug | ✓ |
| C797S | C797 (covalent anchor) | Yes — covalent bonding site | ✓ (via compound profiles) |
| S768I | S768 (αC-helix/activation loop junction) | No — ~12 Å from ligand, allosteric | ✗ |
| L861Q | L861 (activation loop tip) | No — ~14 Å from ligand | ✗ |

**Consequence:** L858R, L861Q, and S768I are excluded. The analysis is anchored on the subset of mutations that physically contact the drug.

---

### v1.1 — Extended contact (8 Å criterion)

**Definition:**  
Same as v1.0 but distance threshold expanded to **8 Å**.  
This captures residues in the second shell of the binding pocket — those that shape the pocket but do not make direct drug contact.

**Expected additions over v1.0:**

| Profile | Key residue(s) | Within 8 Å? | Include? |
|---|---|---|---|
| L858R | L858 | ~10–12 Å — still outside 8 Å threshold | ✗ |
| L861Q | L861 | ~14 Å — still outside | ✗ |
| S768I | S768 | ~12 Å — still outside | ✗ |

*Note: If structural measurement confirms any of these fall within 8 Å, they are automatically included. The criterion is applied mechanically, not by judgment.*

**Consequence:** v1.1 likely produces the same profile list as v1.0 for EGFR. Its value is in confirming that no profiles were incorrectly excluded by the strict 5 Å cutoff.

---

### v2.0 — Functional binding-pocket (structural element criterion)

**Definition:**  
Include a mutation profile if the mutated residue is located in one of the following named structural elements of the EGFR kinase domain, as defined by Yun et al. (2007, PNAS) and the EGFR kinase domain crystal structure annotation:

- **P-loop** (phosphate-binding loop, residues 696–704): G719 is here
- **β3-αC loop** (residues 740–752): Ex19del variants are here
- **Hinge region** (residues 793–800): drug H-bond donor/acceptor
- **Gatekeeper position** (T790): single residue, controls back-pocket access
- **Covalent anchor** (C797): thiol group targeted by irreversible TKIs
- **DFG motif** (residues 855–857): opens/closes the activation loop

**Does NOT include:**
- Activation loop tip (L858, L861) — these mutations destabilize the inactive conformation through a mechanism distal from the drug-binding cleft
- αC-helix (S768) — allosteric, not pocket-lining

**Consequence:** Same profile list as v1.0/v1.1 for current EGFR profiles, but with a named structural element justification rather than a distance cutoff. More readable for a Methods section.

---

### v3.0 — Broad kinase domain (all sensitizing mutations)

**Definition:**  
Include all EGFR kinase domain mutations classified as "sensitizing" or "likely sensitizing" by OncoKB (oncokb.org) regardless of structural location, EXCEPT:  
- Mutations classified as "resistance" only (e.g., T790M alone, C797S alone)  
- Mutations with unknown functional consequence

**Expected additions over v2.0:**

| Profile | OncoKB classification | Include? |
|---|---|---|
| L858R | Sensitizing | ✓ |
| L861Q | Sensitizing | ✓ |
| S768I | Sensitizing (uncommon) | ✓ |

**Consequence:** Widest scope. Includes L858R, L861Q, S768I. Most profiles. This is the "all sensitizing mutations" analysis. Least mechanistically specific — results here are harder to interpret biologically if the correlation weakens.

---

## Decision path

```
Run v1.0
    → p < 0.05 and direction negative?
        YES → publish as pre-specified binding-pocket study (confirmatory framing)
        NO  → run v1.1 (confirm v1.0 profile list is correct)
                → same result → run v2.0
                    → p < 0.05 → publish with v2.0 framing
                    → still no → run v3.0
                        → p < 0.05 → publish with broad sensitizing-mutation framing
                        → still no → update v1_exploratory manuscript
                                      with all findings; publish as exploratory
```

At every step the result — positive or negative — is reported honestly.  
We do not cherry-pick the scope version that gives the best p-value and present it as if it were the only analysis run. All versions tested are reported in a supplementary table.

---

## What to measure for each scope version

For each qualifying profile × drug pair:
1. ΔΔG from GNINA docking (already computed for most pairs; new docking for any missing pairs)
2. KM median TTNT from MSK-CHORD or GENIE BPC
3. n patients and n events

Then: one Spearman correlation per scope version.  
Report: ρ, p, bootstrap 95% CI, n pairs, list of included profiles.
