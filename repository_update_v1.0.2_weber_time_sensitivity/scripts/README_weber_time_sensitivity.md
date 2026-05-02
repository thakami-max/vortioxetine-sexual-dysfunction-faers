# FAERS/AEMS V2 Weber-effect / launch-period time-window sensitivity package

This add-on analysis uses the processed focused V2 primary output to assess whether the main vortioxetine-versus-escitalopram result is driven by early post-launch reporting dynamics.

## Primary time-window sensitivity

The primary time-restricted analysis retains reports with `FDA_DT >= 2017-01-01`, thereby excluding 2014 Q1 through 2016 Q4, the first three full post-approval calendar years for vortioxetine. This follows the proposed Weber-effect sensitivity design while retaining all other primary V2 restrictions:

- Primary suspect (PS) exposure role
- Depression-related indication restriction
- Monotherapy-like restriction
- Narrow sexual-dysfunction preferred-term dictionary
- Active-comparator pairwise case/non-case design

## Optional temporal robustness

The script also computes the primary vortioxetine-versus-escitalopram ROR across three broad calendar periods:

- 2014 Q1-2017 Q4
- 2018 Q1-2021 Q4
- 2022 Q1-2025 Q4

## Installation

Place `faers_v2_weber_time_sensitivity.py` in the main `faers_antidepressant_sd_pipeline` project folder, where the folder `outputs_v2_primary_PS_depression_mono_narrow` is already present.

## Run command

```bash
python3 faers_v2_weber_time_sensitivity.py \
  --v2-out-dir outputs_v2_primary_PS_depression_mono_narrow \
  --out-dir outputs_v2_weber_time_sensitivity \
  --index-drug vortioxetine \
  --comparators escitalopram paroxetine sertraline duloxetine venlafaxine bupropion mirtazapine \
  --steady-start 2017-01-01 \
  --min-cases 3
```

## Key outputs

- `outputs_v2_weber_time_sensitivity/time_window_sensitivity_summary.md`
- `outputs_v2_weber_time_sensitivity/results/time_restricted_pairwise_2017Q1_onward.csv`
- `outputs_v2_weber_time_sensitivity/results/time_restricted_pooled_class_2017Q1_onward.csv`
- `outputs_v2_weber_time_sensitivity/results/time_restricted_vtx_vs_escitalopram_by_sex.csv`
- `outputs_v2_weber_time_sensitivity/results/calendar_period_vtx_vs_escitalopram.csv`
- `outputs_v2_weber_time_sensitivity/results/calendar_period_exposure_counts.csv`

## Interpretation

For pairwise RORs, values below 1 indicate lower sexual-dysfunction reporting with vortioxetine relative to the comparator; values above 1 indicate higher reporting.
