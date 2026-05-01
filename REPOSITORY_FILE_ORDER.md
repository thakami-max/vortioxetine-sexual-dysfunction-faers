# Recommended repository file order

```text
vortioxetine-sexual-dysfunction-faers/
├── README.md
├── requirements.txt
├── CITATION.cff
├── data_download_manifest.csv
├── REPOSITORY_FILE_ORDER.md
├── scripts/
│   ├── faers_v2_active_comparator_pipeline.py
│   ├── faers_v2_print_key_results.py
│   ├── faers_v2_confounder_sensitivity.py
│   └── faers_v2_print_confounder_results.py
├── config/
│   ├── antidepressants_selected_comparators.csv
│   ├── indications_depression.csv
│   ├── sexual_dysfunction_terms_narrow_primary.csv
│   ├── sexual_dysfunction_terms_broad.csv
│   ├── concomitant_drug_groups.csv
│   └── concomitant_condition_terms.csv
├── results_summary/
│   ├── primary_pairwise_2x2_tables.csv
│   ├── prr_chi_square_robustness_table.csv
│   ├── concomitant_flag_counts_by_drug.csv
│   ├── vortioxetine_vs_escitalopram_confounder_exclusion_sensitivities.csv
│   ├── pairwise_after_excluding_core_confounders_except_pde5.csv
│   ├── vortioxetine_vs_escitalopram_stratified_by_core_confounder.csv
│   └── confounder_sensitivity_summary.md
├── figures/
│   ├── figure1_flow_diagram.png
│   └── figure2_forest_plot.png
├── docs/
│   ├── READUS_PV_checklist.docx
│   ├── analysis_commands_and_reproducibility_log.txt
│   ├── README_V2_original.md
│   └── example_commands_V2.md
└── supplementary_material/
    ├── README_ordered_supplementary_material.md
    ├── Supplementary_File_S1_Selected_comparator_drug_dictionary.csv
    ├── Supplementary_File_S2_Depression_indication_dictionary.csv
    ├── Supplementary_File_S3_Narrow_sexual_dysfunction_PT_dictionary.csv
    ├── Supplementary_File_S4_Broad_sexual_dysfunction_PT_dictionary.csv
    ├── Supplementary_File_S5_Complete_primary_2x2_pairwise_case_noncase_tables.csv
    ├── Supplementary_File_S6_PRR_and_chi_square_robustness_table.csv
    ├── Supplementary_File_S7_READUS_PV_checklist.docx
    ├── Supplementary_File_S8_Analysis_commands_and_reproducibility_log/
    ├── Supplementary_File_S9_Contextual_all_antidepressant_analysis/
    ├── Supplementary_File_S10_Concomitant_drug_and_condition_flag_dictionaries/
    ├── Supplementary_File_S11_Concomitant_flag_counts_by_drug.csv
    ├── Supplementary_File_S12_Confounder_exclusion_2x2_tables/
    ├── Supplementary_File_S13_Confounder_stratified_primary_comparison/
    └── Supplementary_File_S14_Figure_files/
```

## Notes

- Do not upload the raw FDA FAERS/AEMS quarterly ZIP files to the repository.
- Instead, use `data_download_manifest.csv` to document the exact public files, source URL, download date, and analytic cut-off.
- Replace placeholder GitHub and Zenodo fields in `README.md` and `CITATION.cff` before public release.
