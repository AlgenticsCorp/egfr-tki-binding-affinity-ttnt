#!/usr/bin/env python3
"""
Patient-level Cox proportional-hazards analysis for Study 2 (EGFR).
Correlates computational binding affinity (ΔΔG) with time-to-next-treatment (TTNT).
Addresses ecological fallacy criticism and performs line-of-therapy sensitivity analysis.

Number-critical analysis: all values derived from local data files, saved to CSV and JSON.

Usage:
    python3 run_patientlevel_cox.py

Output:
    - v2_prespecified/results/patientlevel_dataset.csv (analytical dataset)
    - v2_prespecified/results/patientlevel_cox_results.json (model results)
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

# ============================================================================
# CONFIGURATION
# ============================================================================

# Absolute paths (for re-runability)
REPO_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = REPO_ROOT / 'v1_exploratory' / 'data'
CLINICAL_DATA = REPO_ROOT.parent / 'data' / 'Datasets' / 'msk_chord_2024_clinical_data.tsv'
OUTPUT_DIR = REPO_ROOT / 'v2_prespecified' / 'results'

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Pocket-scope EGFR TKI binding profiles (from preregistration)
POCKET_SCOPE_PROFILES = {
    'Ex19del',
    'G719X',
    'G719X+S768I',
    'G719X+L861Q',
    'L858R+T790M',
    'Ex19del+T790M'
}

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("=" * 80)
print("STEP 1: LOADING DATA")
print("=" * 80)

pfs = pd.read_csv(DATA_DIR / 'pfs' / 'pfs_by_line.csv')
ddg = pd.read_csv(DATA_DIR / 'correlation' / 'merged_ddg_ttnt.csv')
clinical = pd.read_csv(CLINICAL_DATA, sep='\t')

print(f"\nLoaded data:")
print(f"  pfs_by_line.csv: {pfs.shape[0]} rows × {pfs.shape[1]} columns")
print(f"  merged_ddg_ttnt.csv: {ddg.shape[0]} rows × {ddg.shape[1]} columns")
print(f"  msk_chord_2024_clinical_data.tsv: {clinical.shape[0]} rows × {clinical.shape[1]} columns")

# ============================================================================
# STEP 2: BUILD PATIENT-LEVEL DATASET
# ============================================================================

print("\n" + "=" * 80)
print("STEP 2: BUILDING PATIENT-LEVEL DATASET")
print("=" * 80)

# Filter for OSIMERTINIB only
osim = pfs[pfs['drug'].str.upper() == 'OSIMERTINIB'].copy()
print(f"\nOsimertinib lines: {len(osim)} (from {osim['patient_id'].nunique()} unique patients)")

# Filter for pocket-scope profiles
osim_pocket = osim[osim['profile'].isin(POCKET_SCOPE_PROFILES)].copy()
print(f"Pocket-scope osimertinib lines: {len(osim_pocket)} (from {osim_pocket['patient_id'].nunique()} unique patients)")

# Rename columns and code outcomes
osim_pocket = osim_pocket.rename(columns={'ttnt_months': 'duration'})

# Invert censoring flag: in pfs_by_line, censored=0 means event occurred
# Cox model expects event=1 for event, event=0 for censored
osim_pocket['event'] = 1 - osim_pocket['censored']

print(f"\nOutcome coding:")
print(f"  event=1 (progression/next-line): {osim_pocket['event'].sum()}")
print(f"  event=0 (censored): {(1-osim_pocket['event']).sum()}")

# Join ΔΔG by profile+drug
osim_pocket['drug_lower'] = osim_pocket['drug'].str.lower()
ddg['drug_lower'] = ddg['drug'].str.lower()

dataset = osim_pocket.merge(
    ddg[['profile', 'drug_lower', 'ddg']],
    on=['profile', 'drug_lower'],
    how='left'
)

assert dataset['ddg'].notna().all(), "Missing ΔΔG values in merged dataset"

# Get patient-level covariates from clinical data
clinical_dedup = clinical.drop_duplicates(subset=['Patient ID']).copy()
clinical_dedup = clinical_dedup[[
    'Patient ID', 'Current Age', 'Sex', 'Stage (Highest Recorded)'
]].copy()
clinical_dedup = clinical_dedup.rename(columns={
    'Patient ID': 'patient_id',
    'Current Age': 'age',
    'Sex': 'sex',
    'Stage (Highest Recorded)': 'stage'
})

# Join clinical covariates
dataset = dataset.merge(clinical_dedup, on='patient_id', how='left')

# Remove rows with missing covariates
rows_before = len(dataset)
dataset = dataset.dropna(subset=['age', 'sex', 'stage'])
rows_removed = rows_before - len(dataset)

print(f"\nFinal analytical dataset: {len(dataset)} observations from {dataset['patient_id'].nunique()} unique patients")
print(f"  Events: {dataset['event'].sum()}")

# Save dataset
output_dataset_path = OUTPUT_DIR / 'patientlevel_dataset.csv'
dataset.to_csv(output_dataset_path, index=False)

# ============================================================================
# STEP 3: PREPARE FOR MODELING
# ============================================================================

# Create dummy variables for categorical covariates
dataset['sex_male'] = (dataset['sex'] == 'Male').astype(int)
dataset['stage_4'] = (dataset['stage'] == 'Stage 4').astype(int)
dataset['line_group'] = (dataset['line_num'] > 1).astype(int)

# ============================================================================
# STEP 4: FIT COX MODELS
# ============================================================================

print("\n" + "=" * 80)
print("STEP 4: COX PROPORTIONAL-HAZARDS MODELS")
print("=" * 80)

results = {}

# Model 1: Unadjusted
print("\nModel 1: UNADJUSTED (ddg only, clustered by patient)")
cph1 = CoxPHFitter(penalizer=0)
cph1.fit(
    df=dataset,
    duration_col='duration',
    event_col='event',
    formula='ddg',
    cluster_col='patient_id'
)

model1_ddg_coef = cph1.summary.loc['ddg', 'coef']
model1_ddg_hr = np.exp(model1_ddg_coef)
model1_ddg_ci_low = np.exp(cph1.summary.loc['ddg', 'coef lower 95%'])
model1_ddg_ci_high = np.exp(cph1.summary.loc['ddg', 'coef upper 95%'])
model1_ddg_p = cph1.summary.loc['ddg', 'p']
model1_n_events = int(dataset['event'].sum())

results['model_1_unadjusted'] = {
    'name': 'Unadjusted',
    'formula': 'duration ~ ddg',
    'clustering': 'patient_id',
    'ddg_hr': float(model1_ddg_hr),
    'ddg_ci_lower': float(model1_ddg_ci_low),
    'ddg_ci_upper': float(model1_ddg_ci_high),
    'ddg_p_value': float(model1_ddg_p),
    'n': len(dataset),
    'n_events': int(model1_n_events)
}

# Model 2: Primary adjusted
print("Model 2: PRIMARY ADJUSTED (ddg + line + age + sex + stage, clustered by patient)")
cph2 = CoxPHFitter(penalizer=0)
cph2.fit(
    df=dataset,
    duration_col='duration',
    event_col='event',
    formula='ddg + line_num + age + sex_male + stage_4',
    cluster_col='patient_id'
)

model2_ddg_coef = cph2.summary.loc['ddg', 'coef']
model2_ddg_hr = np.exp(model2_ddg_coef)
model2_ddg_ci_low = np.exp(cph2.summary.loc['ddg', 'coef lower 95%'])
model2_ddg_ci_high = np.exp(cph2.summary.loc['ddg', 'coef upper 95%'])
model2_ddg_p = cph2.summary.loc['ddg', 'p']
model2_n_events = int(dataset['event'].sum())

# Schoenfeld proportional-hazards test for the ddg term (primary model)
try:
    from lifelines.statistics import proportional_hazard_test
    ph = proportional_hazard_test(cph2, dataset, time_transform='rank')
    model2_ph_ddg_p = float(ph.summary.loc['ddg', 'p'])
except Exception as e:
    model2_ph_ddg_p = None
    print(f"  (Schoenfeld PH test could not be computed: {e})")

pct_change = (model2_ddg_hr - model1_ddg_hr) / model1_ddg_hr * 100

results['model_2_primary_adjusted'] = {
    'name': 'Primary Adjusted',
    'formula': 'duration ~ ddg + line_num + age + sex_male + stage_4',
    'clustering': 'patient_id',
    'ddg_hr': float(model2_ddg_hr),
    'ddg_ci_lower': float(model2_ddg_ci_low),
    'ddg_ci_upper': float(model2_ddg_ci_high),
    'ddg_p_value': float(model2_ddg_p),
    'n': len(dataset),
    'n_events': int(model2_n_events),
    'schoenfeld_ph_test_p_ddg': model2_ph_ddg_p,
    'percent_change_vs_unadjusted': float(pct_change),
    'covariates_used': ['line_num', 'age', 'sex_male', 'stage_4']
}

# Model 3: Line sensitivity A (same as Model 2)
results['model_3_line_sensitivity_a'] = {
    'name': 'Line Sensitivity A (Adjusted with line)',
    'formula': 'duration ~ ddg + line_num + age + sex_male + stage_4',
    'clustering': 'patient_id',
    'note': 'Same as Model 2 - primary adjusted model',
    'ddg_hr': float(model2_ddg_hr),
    'ddg_ci_lower': float(model2_ddg_ci_low),
    'ddg_ci_upper': float(model2_ddg_ci_high),
    'ddg_p_value': float(model2_ddg_p),
    'n': len(dataset),
    'n_events': int(model2_n_events)
}

# Model 4: Line sensitivity B (stratified by line group)
print("Model 4: LINE SENSITIVITY B (stratified by line group, clustered by patient)")
cph4 = CoxPHFitter(penalizer=0)
cph4.fit(
    df=dataset,
    duration_col='duration',
    event_col='event',
    formula='ddg',
    strata=['line_group'],
    cluster_col='patient_id'
)

model4_ddg_coef = cph4.summary.loc['ddg', 'coef']
model4_ddg_hr = np.exp(model4_ddg_coef)
model4_ddg_ci_low = np.exp(cph4.summary.loc['ddg', 'coef lower 95%'])
model4_ddg_ci_high = np.exp(cph4.summary.loc['ddg', 'coef upper 95%'])
model4_ddg_p = cph4.summary.loc['ddg', 'p']
model4_n_events = int(dataset['event'].sum())

results['model_4_line_sensitivity_b'] = {
    'name': 'Line Sensitivity B (Stratified by line group)',
    'formula': 'duration ~ ddg (stratified by line group)',
    'clustering': 'patient_id',
    'ddg_hr': float(model4_ddg_hr),
    'ddg_ci_lower': float(model4_ddg_ci_low),
    'ddg_ci_upper': float(model4_ddg_ci_high),
    'ddg_p_value': float(model4_ddg_p),
    'n': len(dataset),
    'n_events': int(model4_n_events),
    'strata': 'line_group (frontline vs later)'
}

# ============================================================================
# INTERPRETATION
# ============================================================================

# Direction logic: the group-level correlation is NEGATIVE (ρ=-0.673), i.e. higher ΔΔG
# (weaker binding) → shorter TTNT. The CONFIRMING patient-level result is therefore HR > 1
# per +1 kcal/mol ΔΔG (higher ΔΔG raises the hazard = shortens TTNT). HR < 1 would OPPOSE it.
if model2_ddg_hr > 1:
    direction_agreement = ("AGREE (confirming): HR>1 per +1 kcal/mol ΔΔG means weaker binding "
                           "shortens TTNT — equivalently stronger (more negative) binding lengthens "
                           "TTNT — the same direction as the group-level ρ=-0.673.")
else:
    direction_agreement = ("DISAGREE: HR<1 would imply weaker binding lengthens TTNT, opposite to "
                           "the group-level negative correlation.")

results['interpretation'] = {
    'group_level_correlation': 'ρ = -0.673 (more negative ΔΔG → longer TTNT)',
    'confirming_direction_is': 'HR > 1 per +1 kcal/mol ΔΔG (because the group-level association is negative)',
    'observed_patient_level_primary': f'HR = {model2_ddg_hr:.4f} (p={model2_ddg_p:.4f})',
    'direction_agreement': direction_agreement
}

# ============================================================================
# SAVE AND REPORT
# ============================================================================

print("\n" + "=" * 80)
print("STEP 5: SAVING RESULTS")
print("=" * 80)

output_json = OUTPUT_DIR / 'patientlevel_cox_results.json'
with open(output_json, 'w') as f:
    json.dump(results, f, indent=2)

# Final summary
print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(f"\n{'Model':<45} {'HR (95% CI)':<35} {'p-value':<10} {'N':>6}")
print(f"{'='*100}")
print(f"{'Model 1: Unadjusted':<45} {f'{model1_ddg_hr:.4f} [{model1_ddg_ci_low:.4f},{model1_ddg_ci_high:.4f}]':<35} {model1_ddg_p:>9.4f} {len(dataset):>6}")
print(f"{'Model 2: Primary Adjusted':<45} {f'{model2_ddg_hr:.4f} [{model2_ddg_ci_low:.4f},{model2_ddg_ci_high:.4f}]':<35} {model2_ddg_p:>9.4f} {len(dataset):>6}")
print(f"{'Model 3: Line Sensitivity A':<45} {f'{model2_ddg_hr:.4f} [{model2_ddg_ci_low:.4f},{model2_ddg_ci_high:.4f}]':<35} {model2_ddg_p:>9.4f} {len(dataset):>6}")
print(f"{'Model 4: Line Sensitivity B (Stratified)':<45} {f'{model4_ddg_hr:.4f} [{model4_ddg_ci_low:.4f},{model4_ddg_ci_high:.4f}]':<35} {model4_ddg_p:>9.4f} {len(dataset):>6}")

print(f"\nDataset: {output_dataset_path}")
print(f"Results: {output_json}")
print("\nAnalysis complete.")

