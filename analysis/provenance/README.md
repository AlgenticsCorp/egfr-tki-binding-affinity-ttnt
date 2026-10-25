# Provenance scripts

These two scripts are published as the record of how two reported results were
produced. **They are not runnable from a clean checkout of this repository.**
Their paths refer to the authors' working tree, and both require inputs that are
not redistributable under the source data licences. They are included so that the
procedure behind each number is inspectable, not so that it can be re-executed
without the underlying data.

## `redock_erlotinib_1M17.py`

Re-docks erlotinib against the active-conformation crystal **1M17** (DFG-in/aC-in,
erlotinib-bound AQ4, 2.6 A), replacing the originally-used inactive-conformation
reference 4HJO. This changed the erlotinib wildtype dG from -6.95 to
**-5.39 kcal/mol** and is the reason every erlotinib ddG in `data/` differs from
the pre-2026-07-14 version of this repository.

Requires: a Neurosnap account and `NEUROSNAP_API_KEY` in the environment, the
prepared erlotinib ligand SDF, and the structure-preparation module from the
authors' pipeline. The run was performed once; Neurosnap's hosted GNINA does not
expose seed or exhaustiveness controls, so the erlotinib wildtype dG is a
single-run estimate rather than a 3-seed mean. This methods deviation is
disclosed in the manuscript.

## `run_patientlevel_cox.py`

Fits the patient-level Cox proportional-hazards models reported in Section 3.7 of
the manuscript (ddG as predictor of time to next treatment, with line-of-therapy
and clinical covariate adjustment, clustered on patient).

Requires the MSK-CHORD per-line patient table. **That input is not shipped here.**
MSK-CHORD 2024 is CC BY-NC-ND, so a derived patient-level table is not
redistributed from this repository. The table is reconstructible from cBioPortal
study `msk_chord_2024` (no registration required).

The fitted model output **is** shipped, in `results/patientlevel_cox_results.json`:
576 treatment lines, 404 events, unadjusted ddG HR 1.09 (95% CI 0.99-1.21,
p = 0.080), primary adjusted HR 1.16 (95% CI 1.05-1.28, p = 0.0031).

Dependencies beyond the main script: `pandas`, `lifelines`.
