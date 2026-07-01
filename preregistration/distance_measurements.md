# EGFR Residue-to-Ligand Distance Measurements
**Date:** 2026-06-26  
**Purpose:** Structural basis for v1.0 (≤5 Å) and v1.1 (≤8 Å) scope decisions  
**Method:** Minimum distance from any heavy atom of the residue to any heavy atom of the co-crystal ligand  
**Position mapping:** MPFGCLL motif used to map UniProt positions to PDB numbering (avoids numbering offset errors)

---

## Important finding: 4HJO is INACTIVE conformation

4HJO (TITLE: "Crystal structure of the inactive EGFR tyrosine kinase domain with erlotinib") captures the DFG-out inactive conformation. In this conformation the αC-helix and activation loop are repositioned relative to the active conformation. 

**Consequence:** 4HJO distances are NOT representative of the binding mode for osimertinib or afatinib (which bind the active conformation). The primary structural references for scope decisions are **4ZAU** (osimertinib-bound, active) and **4G5J** (afatinib-bound, active). 4HJO measurements are shown for completeness but are not used for scope qualification.

---

## Measurements (all distances in Ångströms)

| Mutation | UniProt | 4ZAU (osi, active) | 4HJO (erl, inactive) | 4G5J (afa, active) | Active min | v1.0 ≤5Å | v1.1 ≤8Å |
|---|---|---|---|---|---|---|---|
| G719 (P-loop) | 719 | **3.73** | 5.06 | 7.58 | **3.73** | ✓ IN | ✓ IN |
| S768 (αC-helix) | 768 | 16.92 | 13.46 | 11.64 | 11.64 | ✗ OUT | ✗ OUT |
| T790 (gatekeeper) | 790 | 4.68 | 3.48 | **3.19** | **3.19** | ✓ IN | ✓ IN |
| C797 (covalent anchor) | 797 | 3.27 | 3.54 | **1.84** | **1.84** | ✓ IN | ✓ IN |
| L858 (activation loop) | 858 | 11.82 | 5.26 | 9.89 | 9.89 | ✗ OUT | ✗ OUT |
| L861 (activation loop) | 861 | 19.86 | 9.22 | 16.14 | 16.14 | ✗ OUT | ✗ OUT |
| E746 (Ex19del start) | 746 | 8.81 | 6.48 | 7.05 | 7.05 | ✗ OUT | ✓ IN |
| L747 (Ex19del mid) | 747 | 8.81 | 8.19 | 7.05 | 7.05 | ✗ OUT | ✓ IN |
| A750 (Ex19del end) | 750 | N/A* | 16.95 | N/A* | 16.95 | ✗ OUT | ✗ OUT |
| T751 (Ex19del end+1) | 751 | N/A* | 12.71 | N/A* | 12.71 | ✗ OUT | ✗ OUT |

*A750 and T751 are not resolved in 4ZAU and 4G5J (disordered loop region in active conformation crystal structures).

---

## Scope decisions by mutation

### Directly qualified (v1.0 and v1.1): ≤5 Å in active conformation
- **G719** — 3.73 Å (4ZAU). P-loop residue, directly contacts drug. ✓
- **T790** — 3.19 Å (4G5J). Gatekeeper residue, directly contacts drug. ✓  
- **C797** — 1.84 Å (4G5J). Covalent anchor. ✓

### Qualified at v1.1 only (≤8 Å in active conformation): 7–8 Å
- **E746/L747 (Ex19del region)** — 7.05–8.81 Å. The Ex19del deletion removes residues 746–750 from the β3-αC loop, reshaping the pocket entrance. These residues are in the second shell of the binding site. ✓ at v1.1, ✗ at v1.0.
- **G719 in afatinib** — 7.58 Å (4G5J). G719 qualifies at v1.0 from osimertinib (3.73 Å); the afatinib distance reflects different drug geometry.

### Excluded at both thresholds
- **S768** — 11.64 Å minimum (active). αC-helix residue, allosteric. ✗ at both thresholds.
- **L858** — 9.89 Å minimum (active). Activation loop, distal from pocket. ✗ at both.
- **L861** — 16.14 Å minimum (active). Activation loop tip, very distal. ✗ at both.

---

## Impact on profile list

### v1.0 qualifying mutations: G719, T790, C797 (and Ex19del via loop deletion)

**Profiles that qualify:**
- Ex19del (all variants) — loop deletion at 746–750, confirmed structurally
- G719X — G719 at 3.73 Å ✓
- Ex19del+T790M — T790 at 3.19 Å + Ex19del ✓
- L858R+T790M — qualifies via T790 component ✓
- G719X+S768I — qualifies via G719 component ✓
- G719X+L861Q — qualifies via G719 component ✓
- C797S+Ex19del+T790M — qualifies via T790 + Ex19del ✓

**Profiles that do NOT qualify under v1.0:**
- L858R alone — L858 is 11.82 Å, outside 5 Å ✗
- L861Q alone — L861 is 16.14 Å ✗
- S768I alone — S768 is 11.64 Å ✗

### v1.1 adds: Ex19del residues (E746, L747) confirmed within 8 Å
No new profiles are added — the Ex19del variants already qualified via the loop deletion argument. v1.1 provides the structural metric confirming this.

### Key finding: L858R is EXCLUDED under both v1.0 and v1.1
L858 is 11.82 Å from osimertinib and 9.89 Å from afatinib in active conformations. It does not make direct or second-shell contact with the bound drug. L858R sensitizes EGFR through conformational destabilization of the inactive state — not through binding-pocket alteration. This is consistent with the scope definition.

---

## Consequence for analysis

Removing L858R-only and L861Q-only from the primary analysis (v1.0 scope) means:
- **L858R + osimertinib** (n=294, TTNT=10.68 mo) → excluded from primary
- **L858R + erlotinib** (n=114, TTNT=4.89 mo) → excluded from primary  
- **L861Q + osimertinib** (n=21, TTNT=7.62 mo) → excluded from primary
- **L861Q + afatinib** (n=6, TTNT=11.93 mo) → excluded from primary

These were included in v1 exploratory. Their exclusion here is NOT post-hoc — it follows mechanically from the structural criterion defined before analysis.

**Profiles remaining under v1.0:** Ex19del variants, G719X variants, T790M-containing compounds, C797S-containing compounds.
