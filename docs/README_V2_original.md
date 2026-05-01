# FAERS/AEMS V2 focused active-comparator re-analysis package

This package implements the focused V2 re-analysis plan for a clinically interpretable paper:

**Index drug:** vortioxetine  
**Main comparator:** escitalopram  
**Key comparators:** paroxetine, sertraline, duloxetine, venlafaxine  
**Lower-sexual-burden reference drugs:** bupropion, mirtazapine  
**Primary design:** active-comparator, indication-restricted, sex-stratified FAERS/AEMS case/non-case analysis.

The V2 pipeline is designed to be copied into the existing `faers_antidepressant_sd_pipeline/` folder that already contains:

```text
data_raw/
  faers_ascii_2014q1.zip
  ...
  faers_ascii_2025Q4.zip
requirements.txt
```

## Files in this package

```text
faers_v2_active_comparator_pipeline.py
faers_v2_print_key_results.py
config/
  antidepressants_selected_comparators.csv
  sexual_dysfunction_terms_narrow_primary.csv
  indications_depression.csv
README_V2.md
example_commands_V2.md
```

## What is new in V2 compared with the earlier pipeline?

V2 adds:

1. **Selected comparator framework**
   - vortioxetine vs escitalopram
   - vortioxetine vs paroxetine
   - vortioxetine vs sertraline
   - vortioxetine vs duloxetine
   - vortioxetine vs venlafaxine
   - vortioxetine vs bupropion
   - vortioxetine vs mirtazapine
   - pooled SSRI and SNRI context

2. **INDI table processing**
   - Links indications by `PRIMARYID + DRUG_SEQ`.
   - Primary analysis retains only depression-related indications.

3. **Monotherapy-like restriction**
   - Excludes reports with more than one selected antidepressant coded `PS` or `SS`.

4. **Pairwise active-comparator RORs**
   - Pairwise case/non-case RORs are calculated directly.
   - ROR > 1 means higher sexual-dysfunction reporting with vortioxetine versus comparator.
   - ROR < 1 means lower sexual-dysfunction reporting with vortioxetine versus comparator.

5. **Sex-stratified pairwise analyses**
   - Male and female analyses are produced separately.

6. **Optional source/country strata**
   - U.S.-only
   - non-U.S./unknown
   - healthcare-professional-only
   - consumer/non-HCP

## Primary recommended command

Run this from inside your existing `faers_antidepressant_sd_pipeline` folder:

```bash
mkdir -p logs

caffeinate -dimsu python3 faers_v2_active_comparator_pipeline.py \
  --raw-dir data_raw \
  --out-dir outputs_v2_primary_PS_depression_mono_narrow \
  --drug-config config/antidepressants_selected_comparators.csv \
  --event-config config/sexual_dysfunction_terms_narrow_primary.csv \
  --indication-config config/indications_depression.csv \
  --index-drug vortioxetine \
  --comparators escitalopram paroxetine sertraline duloxetine venlafaxine bupropion mirtazapine \
  --role-codes PS \
  --exclude-other-selected-roles PS SS \
  --indication-mode depression \
  --monotherapy-like \
  --run-sex-strata \
  --run-source-strata \
  --min-cases 3 2>&1 | tee logs/v2_primary_PS_depression_mono_narrow.log
```

## Print key results

After the run finishes:

```bash
python3 faers_v2_print_key_results.py --out-dir outputs_v2_primary_PS_depression_mono_narrow
```

## Key output files

```text
outputs_v2_primary_PS_depression_mono_narrow/
  pipeline_v2_summary.md
  processed/
    demo_deduplicated.csv
    selected_antidepressant_exposures_all_scanned_roles.csv
    selected_antidepressant_exposures_primary_roles_pre_indication.csv
    depression_indication_matches.csv
    selected_antidepressant_exposures_final.csv
    sexual_dysfunction_event_matches.csv
    case_level_flags_v2.csv
  results/
    pairwise_drug_comparisons_overall.csv
    pairwise_drug_comparisons_stratified.csv
    pairwise_pooled_class_comparisons_overall.csv
    pairwise_domain_specific_comparisons.csv
    pairwise_pt_specific_comparisons.csv
    selected_exposure_counts_final.csv
    selected_exposure_sd_counts_final.csv
    event_counts.csv
  figures/
    forest_pairwise_vortioxetine_vs_comparators.png
```

## Sensitivity analyses

See `example_commands_V2.md`.

## Interpretation note

This is a spontaneous-reporting case/non-case analysis. RORs are disproportionality metrics and must not be interpreted as incidence, prevalence, relative risk, or proof of causality.
