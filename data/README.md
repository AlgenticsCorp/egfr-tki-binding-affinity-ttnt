# Data

## What is included in this repository

| File | Description | Source |
|---|---|---|
| `merged_ddg_ttnt.csv` | Profile-level ΔΔG (from docking) and median TTNT (from MSK-CHORD KM), n=19 rows | Computed — see analysis script |
| `genie_bpc_profile_summary.csv` | Profile-level median TTNT from GENIE BPC NSCLC v2.0-public, n=34 rows | Computed from GENIE BPC (see below) |

These are aggregate, profile-level summaries with no patient identifiers. They are sufficient to reproduce all primary and secondary correlation results and all four figures.

## What is NOT included

`genie_bpc_cohort.csv` — the patient-level GENIE BPC cohort file (535 rows) used for
Supplementary Table 4 (treatment-line confounding analysis). This file contains
de-identified sample IDs from GENIE BPC and is excluded from this repository in
accordance with the Synapse data access agreement.

To reproduce Supplementary Table 4:

1. Register at [Synapse](https://www.synapse.org) (free).
2. Accept the GENIE BPC data access agreement.
3. Download **GENIE BPC NSCLC v2.0-public** from accession **syn27056697**.
4. Extract the patient-level regimen file and place it here as `genie_bpc_cohort.csv`.
5. Re-run `python analysis/run_analysis.py` — the Supplementary Table 4 section will
   execute automatically when the file is present.

All other results (primary/secondary correlations, LOO analysis, biochemical validation,
GENIE BPC replication) reproduce without this file.

## Raw clinical data sources

| Dataset | Access | Reference |
|---|---|---|
| MSK-CHORD 2024 | Free, no registration — [cBioPortal](https://www.cbioportal.org/study/summary?id=msk_chord_2024) | Memorial Sloan Kettering Cancer Center |
| GENIE BPC NSCLC v2.0-public | Free registration at [Synapse](https://www.synapse.org), accession syn27056697 | AACR Project GENIE Biopharma Collaborative |

The raw patient-level data from these sources were processed by the analysis script
to produce the aggregate CSV files above. Neither source contains patient names or
direct identifiers.
