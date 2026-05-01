# Vortioxetine sexual dysfunction FAERS/AEMS active-comparator analysis

This repository contains analysis scripts, configuration dictionaries, commands, figure files, and summary outputs for the manuscript:

**Sexual Dysfunction Reports Associated With Vortioxetine and Other Antidepressants: A Comparative FAERS/AEMS Disproportionality Analysis**

## Study overview

This was an active-comparator, depression-indication-restricted, sex-stratified FAERS/AEMS case/non-case disproportionality analysis. The index drug was **vortioxetine**. The prespecified primary comparator was **escitalopram**; secondary comparators were **paroxetine**, **sertraline**, **duloxetine**, **venlafaxine**, **bupropion**, and **mirtazapine**. Contextual class comparisons included pooled SSRIs and pooled SNRIs.

The primary analysis used:

- Primary-suspect (PS) drug role
- Depression-related indication restriction using linked INDI records
- A monotherapy-like restriction excluding reports with more than one selected antidepressant coded PS or secondary suspect (SS)
- A narrow sexual-dysfunction preferred-term dictionary
- Case/non-case pairwise reporting odds ratios (RORs)

## Public data source

The analysis used public FDA FAERS/AEMS quarterly ASCII extract files from **2014 Q1 through 2025 Q4**.

- Data download date: **2026-03-31**
- Prespecified analytic cut-off: **2025 Q4**
- Raw public files: FDA FAERS/AEMS quarterly ASCII extracts
- Raw quarterly ZIP files are **not redistributed** in this repository because they are publicly available from the FDA and are large.

A file-level manifest is provided in `data_download_manifest.csv`.

## Repository contents

```text
README.md
requirements.txt
CITATION.cff
data_download_manifest.csv
REPOSITORY_FILE_ORDER.md
scripts/
config/
results_summary/
figures/
docs/
supplementary_material/
```

## Software environment

The submitted analysis was run using:

- Python 3.9.x
- pandas 2.3.3
- NumPy 2.0.2
- Matplotlib 3.9.4

The Python package requirements are provided in `requirements.txt`.

## Reproducibility steps

1. Download the public FDA FAERS/AEMS ASCII quarterly files listed in `data_download_manifest.csv`.
2. Place the raw ZIP files in a local folder named `data_raw/`.
3. Install the Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

4. Run the primary V2 analysis using the commands in `docs/analysis_commands_and_reproducibility_log.txt`.

The main scripts are located in `scripts/`:

```text
scripts/faers_v2_active_comparator_pipeline.py
scripts/faers_v2_print_key_results.py
scripts/faers_v2_confounder_sensitivity.py
scripts/faers_v2_print_confounder_results.py
```

## Supplementary material

The folder `supplementary_material/` contains the ordered Supplementary Files S1-S14 used for submission:

- S1 selected comparator drug dictionary
- S2 depression-indication dictionary
- S3 narrow sexual-dysfunction preferred-term dictionary
- S4 broad sexual-dysfunction preferred-term dictionary
- S5 complete primary 2 x 2 pairwise case/non-case tables
- S6 PRR and chi-square robustness table
- S7 READUS-PV checklist
- S8 analysis commands, reproducibility log, and scripts
- S9 contextual broader all-antidepressant analysis
- S10 concomitant drug and condition flag dictionaries
- S11 concomitant flag counts by drug
- S12 confounder-exclusion 2 x 2 tables
- S13 stratified analyses for the primary vortioxetine-versus-escitalopram comparison
- S14 figure source files

## Interpretation note

FAERS/AEMS is a spontaneous reporting system. The analyses in this repository are intended for signal detection and clinical-contextual interpretation. They do not estimate incidence, prevalence, comparative risk, or causality.

## Citation

Please cite the associated manuscript and the archived repository release. After creating a Zenodo archive, replace the placeholder below with the final DOI:

```text
Hakami, T., et al. (2026). Vortioxetine sexual dysfunction FAERS/AEMS active-comparator analysis: scripts and configuration files (v1.0.0). Zenodo. DOI: INSERT_DOI
```
