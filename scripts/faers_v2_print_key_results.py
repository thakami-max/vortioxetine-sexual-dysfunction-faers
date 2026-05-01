#!/usr/bin/env python3
"""Print key V2 FAERS/AEMS active-comparator results after a completed run."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def fmt_table(df: pd.DataFrame, cols: list[str]) -> str:
    return df[cols].to_string(index=False)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    out = args.out_dir

    overall_path = out / "results" / "pairwise_drug_comparisons_overall.csv"
    strat_path = out / "results" / "pairwise_drug_comparisons_stratified.csv"
    pooled_path = out / "results" / "pairwise_pooled_class_comparisons_overall.csv"
    domain_path = out / "results" / "pairwise_domain_specific_comparisons.csv"
    counts_path = out / "results" / "selected_exposure_sd_counts_final.csv"

    print("\n=== Summary ===")
    summary = out / "pipeline_v2_summary.md"
    if summary.exists():
        print(summary.read_text())

    if counts_path.exists():
        print("\n=== Final selected exposure counts ===")
        df = pd.read_csv(counts_path)
        print(fmt_table(df, ["drug_class", "generic", "n_exposed", "sd_cases", "sd_reporting_proportion"]))

    if overall_path.exists():
        print("\n=== Pairwise drug comparisons: overall ===")
        df = pd.read_csv(overall_path)
        df = df.sort_values(["comparator"])
        cols = ["index_drug", "comparator", "a_index_with_SD", "n_index", "c_comparator_with_SD", "n_comparator", "ror", "ror025", "ror975", "signal_higher_vortioxetine_lower95_gt1", "signal_lower_vortioxetine_upper95_lt1"]
        print(fmt_table(df, cols))

    if pooled_path.exists():
        print("\n=== Pooled class comparisons ===")
        df = pd.read_csv(pooled_path)
        cols = ["index_drug", "comparator", "a_index_with_SD", "n_index", "c_comparator_with_SD", "n_comparator", "ror", "ror025", "ror975"]
        print(fmt_table(df, cols))

    if strat_path.exists():
        print("\n=== Sex-stratified pairwise comparisons ===")
        df = pd.read_csv(strat_path)
        df = df[df["stratum"].eq("sex_std") & df["stratum_value"].isin(["Male", "Female"])]
        cols = ["stratum_value", "index_drug", "comparator", "a_index_with_SD", "n_index", "c_comparator_with_SD", "n_comparator", "ror", "ror025", "ror975"]
        print(fmt_table(df, cols))

    if domain_path.exists():
        print("\n=== Domain-specific pairwise comparisons: vortioxetine vs escitalopram ===")
        df = pd.read_csv(domain_path)
        df = df[df["comparator"].eq("escitalopram")].sort_values(["a_index_with_event", "ror"], ascending=[False, False])
        cols = ["domain", "index_drug", "comparator", "a_index_with_event", "c_comparator_with_event", "ror", "ror025", "ror975"]
        print(fmt_table(df, cols))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
