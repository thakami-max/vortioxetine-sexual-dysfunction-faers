# Example V2 commands

Run all commands from inside your existing `faers_antidepressant_sd_pipeline` folder after copying these V2 files into it.

## 1. Primary V2 analysis

PS-only, depression-indication restricted, monotherapy-like, narrow sexual dysfunction endpoint.

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

Print results:

```bash
python3 faers_v2_print_key_results.py --out-dir outputs_v2_primary_PS_depression_mono_narrow
```

## 2. PS + SS sensitivity analysis

```bash
caffeinate -dimsu python3 faers_v2_active_comparator_pipeline.py \
  --raw-dir data_raw \
  --out-dir outputs_v2_sensitivity_PS_SS_depression_mono_narrow \
  --drug-config config/antidepressants_selected_comparators.csv \
  --event-config config/sexual_dysfunction_terms_narrow_primary.csv \
  --indication-config config/indications_depression.csv \
  --index-drug vortioxetine \
  --comparators escitalopram paroxetine sertraline duloxetine venlafaxine bupropion mirtazapine \
  --role-codes PS SS \
  --exclude-other-selected-roles PS SS \
  --indication-mode depression \
  --monotherapy-like \
  --run-sex-strata \
  --run-source-strata \
  --min-cases 3 2>&1 | tee logs/v2_sensitivity_PS_SS_depression_mono_narrow.log
```

## 3. All-indications sensitivity analysis

```bash
caffeinate -dimsu python3 faers_v2_active_comparator_pipeline.py \
  --raw-dir data_raw \
  --out-dir outputs_v2_sensitivity_PS_all_indications_mono_narrow \
  --drug-config config/antidepressants_selected_comparators.csv \
  --event-config config/sexual_dysfunction_terms_narrow_primary.csv \
  --index-drug vortioxetine \
  --comparators escitalopram paroxetine sertraline duloxetine venlafaxine bupropion mirtazapine \
  --role-codes PS \
  --exclude-other-selected-roles PS SS \
  --indication-mode none \
  --monotherapy-like \
  --run-sex-strata \
  --run-source-strata \
  --min-cases 3 2>&1 | tee logs/v2_sensitivity_PS_all_indications_mono_narrow.log
```

## 4. Broad sexual dysfunction endpoint sensitivity analysis

```bash
caffeinate -dimsu python3 faers_v2_active_comparator_pipeline.py \
  --raw-dir data_raw \
  --out-dir outputs_v2_sensitivity_PS_depression_mono_broad \
  --drug-config config/antidepressants_selected_comparators.csv \
  --event-config config/sexual_dysfunction_terms_narrow_primary.csv \
  --indication-config config/indications_depression.csv \
  --index-drug vortioxetine \
  --comparators escitalopram paroxetine sertraline duloxetine venlafaxine bupropion mirtazapine \
  --role-codes PS \
  --exclude-other-selected-roles PS SS \
  --indication-mode depression \
  --monotherapy-like \
  --include-sensitivity-terms \
  --run-sex-strata \
  --run-source-strata \
  --min-cases 3 2>&1 | tee logs/v2_sensitivity_PS_depression_mono_broad.log
```
