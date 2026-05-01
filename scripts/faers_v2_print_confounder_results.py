#!/usr/bin/env python3
"""Print key results from the V2 concomitant drug/condition sensitivity analysis."""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def fmt_table(df: pd.DataFrame, cols: list[str]) -> str:
    if df.empty:
        return "(no rows)"
    show = df[cols].copy()
    for c in ["ror", "ror025", "ror975", "prr", "chi2"]:
        if c in show.columns:
            show[c] = show[c].astype(float).map(lambda x: f"{x:.2f}")
    return show.to_string(index=False)


def main() -> int:
    p = argparse.ArgumentParser(description="Print key confounder sensitivity results")
    p.add_argument("--out-dir", default="outputs_v2_confounder_sensitivity")
    p.add_argument("--comparator", default="escitalopram")
    args = p.parse_args()
    out = Path(args.out_dir)
    summary = out / "confounder_sensitivity_summary.md"
    overall_path = out / "results" / "pairwise_confounder_sensitivity_overall.csv"
    flags_path = out / "results" / "confounder_flag_counts_by_drug.csv"
    strat_path = out / "results" / "pairwise_confounder_sensitivity_stratified.csv"

    if summary.exists():
        print("\n=== Summary ===")
        print(summary.read_text(encoding="utf-8"))

    overall = pd.read_csv(overall_path)
    flags = pd.read_csv(flags_path)

    print("\n=== Confounder flag counts by drug (compact) ===")
    flag_cols = [c for c in flags.columns if c.endswith("_cases")]
    compact_cols = ["generic", "n_exposed", "sd_cases"] + flag_cols[:12]
    print(fmt_table(flags.sort_values("generic"), compact_cols))

    print(f"\n=== Vortioxetine vs {args.comparator}: exclusions/sensitivities ===")
    esc = overall[overall["comparator"].eq(args.comparator)].copy()
    order = [
        "primary_v2_replicated",
        "exclude_core_confounders_except_pde5",
        "exclude_any_confounder_except_pde5",
        "exclude_any_confounder_including_pde5",
        "exclude_antipsychotic",
        "exclude_five_ari_alpha_blocker",
        "exclude_opioid",
        "exclude_hormonal_therapy",
        "exclude_urologic_endocrine_condition",
        "exclude_pde5_or_ed_treatment",
    ]
    esc["_ord"] = esc["analysis"].apply(lambda x: order.index(x) if x in order else 999)
    esc = esc.sort_values(["_ord", "analysis"])
    cols = ["analysis", "a_index_with_SD", "n_index", "c_comparator_with_SD", "n_comparator", "ror", "ror025", "ror975", "direction"]
    print(fmt_table(esc, cols))

    print("\n=== All pairwise comparisons after excluding core confounders except PDE5 ===")
    core = overall[overall["analysis"].eq("exclude_core_confounders_except_pde5")].sort_values("comparator")
    cols2 = ["comparator", "a_index_with_SD", "n_index", "c_comparator_with_SD", "n_comparator", "ror", "ror025", "ror975", "direction"]
    print(fmt_table(core, cols2))

    if strat_path.exists():
        strat = pd.read_csv(strat_path)
        print("\n=== Vortioxetine vs escitalopram stratified by core confounder ===")
        s = strat[(strat["comparator"].eq(args.comparator)) & (strat["stratum"].eq("any_core_confounder"))].copy()
        cols3 = ["stratum", "stratum_value", "a_index_with_SD", "n_index", "c_comparator_with_SD", "n_comparator", "ror", "ror025", "ror975", "direction"]
        print(fmt_table(s.sort_values("stratum_value"), cols3))

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
