# Data Attribution and Licensing

This directory contains three data files. Their licensing status differs:

---

## ddg_values.csv — MIT License (free to use, modify, redistribute)

**What it contains:** Computational ΔΔG binding affinity values for EGFR mutation–TKI
pairs, derived from AutoDock Vina molecular docking simulations performed by the
study authors.

**Source:** Original computational work by Ahmad Abboud, Algentics Research.  
**License:** [MIT](../LICENSE)  
**No restrictions.** This file contains only the authors' own computations.

---

## merged_ddg_ttnt.csv — CC BY-NC-ND 4.0 (non-commercial, no adapted derivatives)

**What it contains:** Profile-level Kaplan–Meier median TTNT (months) and patient
counts for each mutation–drug pair, derived from MSK-CHORD 2024. Also contains the
`ddg` column (same values as `ddg_values.csv`).

**Derived from:** MSK-CHORD 2024, Memorial Sloan Kettering Cancer Center  
**Raw data access:** cBioPortal — https://www.cbioportal.org/study/summary?id=msk_chord_2024  
**License of raw data:** [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/)  
**Attribution:** Memorial Sloan Kettering Cancer Center. MSK-CHORD 2024.
cBioPortal for Cancer Genomics, 2024.

The TTNT medians and patient counts in this file are aggregate, profile-level summary
statistics (19 rows). No patient identifiers are present. These values are identical
to Table 1 of the associated manuscript and represent the scientific results of the
analysis, not a redistribution of the original patient-level dataset.

**Non-commercial use only.** For commercial applications, contact
datarequests@mskcc.org.

---

## genie_bpc_profile_summary.csv — CC BY-NC-ND 4.0 (non-commercial, no adapted derivatives)

**What it contains:** Profile-level KM median TTNT (months), patient counts, and
event counts for each mutation–drug pair in the GENIE BPC NSCLC cohort, used for the
independent replication analysis (Figure 4 and secondary correlation).

**Derived from:** AACR Project GENIE Biopharma Collaborative, GENIE BPC NSCLC
v2.0-public (Synapse accession syn27056697)  
**Raw data access:** [Synapse syn27056697](https://www.synapse.org/#!Synapse:syn27056697)  
**License of raw data:** [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/)  
**Attribution:** AACR Project GENIE Biopharma Collaborative. GENIE BPC NSCLC
v2.0-public. Synapse syn27056697, 2022.

The TTNT medians in this file are aggregate, profile-level summary statistics (33 rows).
No patient identifiers are present. The underlying patient-level data require a Synapse
data access agreement (see `README.md` for download instructions).

**Non-commercial use only.** For commercial applications, contact the AACR
(genie@aacr.org).

---

## Summary

| File | Contains PHI? | License | Redistribution |
|---|:---:|---|---|
| `ddg_values.csv` | No | MIT | Unrestricted |
| `merged_ddg_ttnt.csv` | No | CC BY-NC-ND 4.0 | Non-commercial only |
| `genie_bpc_profile_summary.csv` | No | CC BY-NC-ND 4.0 | Non-commercial only |
| `genie_bpc_cohort.csv` | Sample IDs | GENIE BPC DAA | Not included — download from Synapse |

All files are aggregate statistics; none contain patient names, dates of birth,
addresses, or other direct identifiers.
