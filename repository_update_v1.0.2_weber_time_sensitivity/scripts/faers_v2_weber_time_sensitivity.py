#!/usr/bin/env python3
"""
FAERS/AEMS V2 Weber-effect / launch-period time-window sensitivity analysis.

Purpose
-------
Uses the already processed focused V2 primary output to test whether the primary
vortioxetine-versus-escitalopram finding is driven by launch-period reporting dynamics.

Primary time-window sensitivity:
  - Restricts the primary V2 cohort to reports with FDA_DT >= 2017-01-01.
  - This excludes 2014 Q1-2016 Q4, the first three full post-approval calendar years.

Optional temporal robustness:
  - Estimates vortioxetine-vs-escitalopram RORs across 2014-2017, 2018-2021, and 2022-2025.

Inputs are read from:
  <v2-out-dir>/processed/selected_antidepressant_exposures_final.csv
  <v2-out-dir>/processed/case_level_flags_v2.csv

Outputs are written to:
  <out-dir>/results/*.csv
  <out-dir>/time_window_sensitivity_summary.md
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


def parse_date_col(s: pd.Series) -> pd.Series:
    raw = s.astype("string").fillna("").str.replace(r"[^0-9]", "", regex=True)
    dt1 = pd.to_datetime(raw, format="%Y%m%d", errors="coerce")
    dt2 = pd.to_datetime(s, errors="coerce")
    return dt1.fillna(dt2)


def find_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    # fuzzy contains fallback
    for cand in candidates:
        c0 = cand.lower()
        for low, orig in lower.items():
            if c0 in low:
                return orig
    return None


def require_column(df: pd.DataFrame, candidates: Sequence[str], label: str) -> str:
    col = find_column(df, candidates)
    if col is None:
        raise ValueError(f"Could not find {label}; tried {candidates}; available columns: {list(df.columns)}")
    return col


def normalize_generic(x: object) -> str:
    return str(x).strip().lower()


def normalize_class(x: object) -> str:
    return str(x).strip().lower()


def standardize_sex_value(x: object) -> str:
    v = str(x).strip().upper()
    if v in {"M", "MALE"}:
        return "Male"
    if v in {"F", "FEMALE"}:
        return "Female"
    return "Unknown"


def stats_from_cells(a: int, b: int, c: int, d: int) -> Dict[str, object]:
    aa, bb, cc, dd = map(float, [a, b, c, d])
    corrected = False
    if min(aa, bb, cc, dd) == 0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
        corrected = True
    ror = (aa * dd) / (bb * cc) if bb * cc > 0 else np.nan
    se = math.sqrt((1 / aa) + (1 / bb) + (1 / cc) + (1 / dd)) if min(aa, bb, cc, dd) > 0 else np.nan
    lo = math.exp(math.log(ror) - 1.96 * se) if ror > 0 and not np.isnan(se) else np.nan
    hi = math.exp(math.log(ror) + 1.96 * se) if ror > 0 and not np.isnan(se) else np.nan
    prr = (aa / (aa + bb)) / (cc / (cc + dd)) if (aa + bb) > 0 and (cc + dd) > 0 and cc > 0 else np.nan
    n = aa + bb + cc + dd
    denom = (aa + bb) * (cc + dd) * (aa + cc) * (bb + dd)
    chi2 = (n * ((aa * dd - bb * cc) ** 2) / denom) if denom > 0 else np.nan
    if not np.isnan(ror):
        if hi < 1:
            direction = "lower_with_vortioxetine"
        elif lo > 1:
            direction = "higher_with_vortioxetine"
        else:
            direction = "no_clear_difference"
    else:
        direction = "not_estimable"
    return {
        "ror": ror,
        "ror025": lo,
        "ror975": hi,
        "prr": prr,
        "chi2": chi2,
        "haldane_correction": corrected,
        "direction": direction,
    }


def case_sets_by_generic(exposures: pd.DataFrame) -> Dict[str, Set[str]]:
    return {g: set(x["primaryid"].astype(str)) for g, x in exposures.groupby("generic")}


def compute_pairwise(
    exposures: pd.DataFrame,
    outcome_cases: Set[str],
    index_drug: str,
    comparator_label: str,
    comparator_cases: Set[str],
    min_cases: int = 3,
    stratum: str = "overall",
    stratum_value: str = "overall",
) -> Dict[str, object]:
    generic_sets = case_sets_by_generic(exposures)
    index_cases = set(generic_sets.get(index_drug.lower(), set()))
    comp_cases = set(comparator_cases)
    overlap = index_cases & comp_cases
    if overlap:
        index_cases -= overlap
        comp_cases -= overlap
    a = len(index_cases & outcome_cases)
    b = len(index_cases - outcome_cases)
    c = len(comp_cases & outcome_cases)
    d = len(comp_cases - outcome_cases)
    stats = stats_from_cells(a, b, c, d)
    row = {
        "stratum": stratum,
        "stratum_value": stratum_value,
        "index_drug": index_drug.lower(),
        "comparator": comparator_label,
        "a_index_with_SD": a,
        "b_index_without_SD": b,
        "c_comparator_with_SD": c,
        "d_comparator_without_SD": d,
        "n_index": a + b,
        "n_comparator": c + d,
        "sd_reporting_prop_index": (a / (a + b)) if (a + b) else np.nan,
        "sd_reporting_prop_comparator": (c / (c + d)) if (c + d) else np.nan,
        **stats,
        "signal_higher_vortioxetine_lower95_gt1": bool(a >= min_cases and stats["ror025"] > 1),
        "signal_lower_vortioxetine_upper95_lt1": bool(c >= min_cases and stats["ror975"] < 1),
    }
    return row


def compute_drug_pairwise(
    exposures: pd.DataFrame,
    outcome_cases: Set[str],
    index_drug: str,
    comparators: Sequence[str],
    min_cases: int = 3,
    stratum: str = "overall",
    stratum_value: str = "overall",
) -> pd.DataFrame:
    generic_sets = case_sets_by_generic(exposures)
    rows: List[Dict[str, object]] = []
    for comp in comparators:
        rows.append(compute_pairwise(
            exposures, outcome_cases, index_drug, comp.lower(), generic_sets.get(comp.lower(), set()),
            min_cases=min_cases, stratum=stratum, stratum_value=stratum_value
        ))
    return pd.DataFrame(rows)


def compute_pooled_class(
    exposures: pd.DataFrame,
    outcome_cases: Set[str],
    index_drug: str,
    class_labels: Sequence[str] = ("SSRI", "SNRI"),
    min_cases: int = 3,
    stratum: str = "overall",
    stratum_value: str = "overall",
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for cls in class_labels:
        comp_cases = set(exposures.loc[
            (exposures["generic"] != index_drug.lower()) &
            (exposures["drug_class"].str.lower() == cls.lower()),
            "primaryid"
        ].astype(str))
        rows.append(compute_pairwise(
            exposures, outcome_cases, index_drug, f"pooled_{cls.lower()}", comp_cases,
            min_cases=min_cases, stratum=stratum, stratum_value=stratum_value
        ))
    return pd.DataFrame(rows)


def load_inputs(v2_out_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    exposures_path = v2_out_dir / "processed" / "selected_antidepressant_exposures_final.csv"
    flags_path = v2_out_dir / "processed" / "case_level_flags_v2.csv"
    if not exposures_path.exists():
        raise FileNotFoundError(f"Missing {exposures_path}")
    if not flags_path.exists():
        raise FileNotFoundError(f"Missing {flags_path}")

    exposures = pd.read_csv(exposures_path, dtype=str, low_memory=False)
    flags = pd.read_csv(flags_path, dtype=str, low_memory=False)

    primary_col_e = require_column(exposures, ["primaryid", "primary_id"], "primaryid in exposures")
    generic_col = require_column(exposures, ["generic", "drug", "exposure"], "generic drug column in exposures")
    class_col = require_column(exposures, ["drug_class", "class"], "drug class column in exposures")
    exposures = exposures.rename(columns={primary_col_e: "primaryid", generic_col: "generic", class_col: "drug_class"}).copy()
    exposures["primaryid"] = exposures["primaryid"].astype(str)
    exposures["generic"] = exposures["generic"].map(normalize_generic)
    exposures["drug_class"] = exposures["drug_class"].astype(str).str.strip()

    primary_col_f = require_column(flags, ["primaryid", "primary_id"], "primaryid in case flags")
    fda_col = require_column(flags, ["fda_dt", "fda_date", "receiptdate", "receipt_date"], "FDA receipt date in case flags")
    outcome_col = require_column(flags, ["sd_any", "has_sd", "sexual_dysfunction", "sd_case"], "sexual dysfunction flag in case flags")
    sex_col = find_column(flags, ["sex", "sex_std", "sex_standardized"])
    flags = flags.rename(columns={primary_col_f: "primaryid", fda_col: "fda_dt", outcome_col: "sd_any"}).copy()
    flags["primaryid"] = flags["primaryid"].astype(str)
    flags["fda_dt_parsed"] = parse_date_col(flags["fda_dt"])
    flags["sd_any_bool"] = flags["sd_any"].astype(str).str.upper().isin(["TRUE", "1", "YES", "Y", "T"])
    if sex_col and sex_col != "sex":
        flags["sex"] = flags[sex_col]
    elif "sex" not in flags.columns:
        flags["sex"] = "Unknown"
    flags["sex_std"] = flags["sex"].map(standardize_sex_value)

    # Attach date/sex/outcome to exposures. Keep only exposures present in case flags.
    keep = flags[["primaryid", "fda_dt_parsed", "sd_any_bool", "sex_std"]].drop_duplicates("primaryid")
    exposures = exposures.merge(keep, on="primaryid", how="inner")
    return exposures, flags


def subset_by_date(exposures: pd.DataFrame, start: Optional[str], end: Optional[str]) -> pd.DataFrame:
    mask = pd.Series(True, index=exposures.index)
    if start:
        mask &= exposures["fda_dt_parsed"] >= pd.Timestamp(start)
    if end:
        mask &= exposures["fda_dt_parsed"] <= pd.Timestamp(end)
    return exposures[mask].copy()


def outcome_cases_from_exposures(exposures: pd.DataFrame) -> Set[str]:
    return set(exposures.loc[exposures["sd_any_bool"], "primaryid"].astype(str))


def summarize_window(label: str, exposures: pd.DataFrame, index_drug: str, comparators: Sequence[str]) -> Dict[str, object]:
    out = {"window": label, "n_exposure_cases": exposures["primaryid"].nunique()}
    for drug in [index_drug] + list(comparators):
        x = exposures[exposures["generic"] == drug.lower()]
        out[f"n_{drug.lower()}"] = x["primaryid"].nunique()
        out[f"sd_{drug.lower()}"] = x.loc[x["sd_any_bool"], "primaryid"].nunique()
    return out


def write_summary(out_dir: Path, overall_summary: Dict[str, object], primary: pd.DataFrame, pooled: pd.DataFrame, sex: pd.DataFrame, periods: pd.DataFrame) -> None:
    lines: List[str] = []
    lines.append("# Weber-effect / launch-period time-window sensitivity analysis summary")
    lines.append("")
    lines.append("## Design")
    lines.append("The time-window sensitivity analysis used the processed primary V2 cohort and restricted reports by FDA_DT, the FDA receipt date used in the case-versioning/deduplication logic.")
    lines.append("The main time-restricted analysis excluded 2014 Q1-2016 Q4 and retained reports with FDA_DT in 2017 Q1 or later.")
    lines.append("")
    lines.append("## Time-restricted cohort")
    for k, v in overall_summary.items():
        lines.append(f"- {k}: {v:,}" if isinstance(v, int) else f"- {k}: {v}")
    lines.append("")
    lines.append("## Vortioxetine vs escitalopram, 2017 Q1 onward")
    vte = primary[primary["comparator"].eq("escitalopram")]
    if not vte.empty:
        r = vte.iloc[0]
        lines.append(f"- ROR: {r['ror']:.3f} (95% CI {r['ror025']:.3f}-{r['ror975']:.3f})")
        lines.append(f"- Counts: vortioxetine {int(r['a_index_with_SD'])}/{int(r['n_index'])}; escitalopram {int(r['c_comparator_with_SD'])}/{int(r['n_comparator'])}")
        lines.append(f"- Direction: {r['direction']}")
    lines.append("")
    lines.append("## Pairwise results, 2017 Q1 onward")
    cols = ["comparator", "a_index_with_SD", "n_index", "c_comparator_with_SD", "n_comparator", "ror", "ror025", "ror975", "direction"]
    lines.append(primary[cols].to_string(index=False))
    lines.append("")
    lines.append("## Pooled class context, 2017 Q1 onward")
    lines.append(pooled[cols].to_string(index=False))
    lines.append("")
    lines.append("## Sex-stratified vortioxetine vs escitalopram, 2017 Q1 onward")
    lines.append(sex[cols + ["stratum", "stratum_value"]].to_string(index=False))
    lines.append("")
    lines.append("## Calendar-period analysis, vortioxetine vs escitalopram")
    period_cols = ["stratum_value", "a_index_with_SD", "n_index", "c_comparator_with_SD", "n_comparator", "ror", "ror025", "ror975", "direction"]
    lines.append(periods[period_cols].to_string(index=False))
    (out_dir / "time_window_sensitivity_summary.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="V2 FAERS/AEMS time-window sensitivity analysis for launch-period/Weber effect")
    p.add_argument("--v2-out-dir", type=Path, default=Path("outputs_v2_primary_PS_depression_mono_narrow"))
    p.add_argument("--out-dir", type=Path, default=Path("outputs_v2_weber_time_sensitivity"))
    p.add_argument("--index-drug", default="vortioxetine")
    p.add_argument("--comparators", nargs="+", default=["escitalopram", "paroxetine", "sertraline", "duloxetine", "venlafaxine", "bupropion", "mirtazapine"])
    p.add_argument("--steady-start", default="2017-01-01", help="Start date for steady-state sensitivity; default excludes 2014-2016")
    p.add_argument("--min-cases", type=int, default=3)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "results").mkdir(parents=True, exist_ok=True)

    exposures, flags = load_inputs(args.v2_out_dir)
    all_n = exposures["primaryid"].nunique()
    all_sd = exposures.loc[exposures["sd_any_bool"], "primaryid"].nunique()

    steady = subset_by_date(exposures, args.steady_start, None)
    steady_outcome = outcome_cases_from_exposures(steady)

    pairwise = compute_drug_pairwise(steady, steady_outcome, args.index_drug, args.comparators, args.min_cases,
                                     stratum="time_window", stratum_value=f">={args.steady_start}")
    pooled = compute_pooled_class(steady, steady_outcome, args.index_drug, ["SSRI", "SNRI"], args.min_cases,
                                  stratum="time_window", stratum_value=f">={args.steady_start}")

    # Sex-stratified primary comparator only.
    sex_rows: List[pd.DataFrame] = []
    for sexval in ["Male", "Female"]:
        sub = steady[steady["sex_std"] == sexval].copy()
        outcome = outcome_cases_from_exposures(sub)
        sex_rows.append(compute_drug_pairwise(sub, outcome, args.index_drug, ["escitalopram"], args.min_cases,
                                              stratum="sex", stratum_value=sexval))
    sex_df = pd.concat(sex_rows, ignore_index=True) if sex_rows else pd.DataFrame()

    # Three period stratified analysis for primary comparator.
    periods = [
        ("2014Q1-2017Q4", "2014-01-01", "2017-12-31"),
        ("2018Q1-2021Q4", "2018-01-01", "2021-12-31"),
        ("2022Q1-2025Q4", "2022-01-01", "2025-12-31"),
    ]
    period_rows: List[pd.DataFrame] = []
    period_summary_rows: List[Dict[str, object]] = []
    for label, start, end in periods:
        sub = subset_by_date(exposures, start, end)
        outcome = outcome_cases_from_exposures(sub)
        res = compute_drug_pairwise(sub, outcome, args.index_drug, ["escitalopram"], args.min_cases,
                                    stratum="calendar_period", stratum_value=label)
        period_rows.append(res)
        period_summary_rows.append(summarize_window(label, sub, args.index_drug, args.comparators))
    period_df = pd.concat(period_rows, ignore_index=True) if period_rows else pd.DataFrame()
    period_summary = pd.DataFrame(period_summary_rows)

    # Counts for main steady subset.
    overall_summary = summarize_window(f">={args.steady_start}", steady, args.index_drug, args.comparators)
    overall_summary["all_primary_v2_exposure_cases"] = all_n
    overall_summary["all_primary_v2_sd_cases"] = all_sd

    # Save outputs.
    pairwise.to_csv(args.out_dir / "results" / "time_restricted_pairwise_2017Q1_onward.csv", index=False)
    pooled.to_csv(args.out_dir / "results" / "time_restricted_pooled_class_2017Q1_onward.csv", index=False)
    sex_df.to_csv(args.out_dir / "results" / "time_restricted_vtx_vs_escitalopram_by_sex.csv", index=False)
    period_df.to_csv(args.out_dir / "results" / "calendar_period_vtx_vs_escitalopram.csv", index=False)
    period_summary.to_csv(args.out_dir / "results" / "calendar_period_exposure_counts.csv", index=False)
    pd.DataFrame([overall_summary]).to_csv(args.out_dir / "results" / "time_restricted_exposure_counts.csv", index=False)

    write_summary(args.out_dir, overall_summary, pairwise, pooled, sex_df, period_df)

    print(f"Done. Outputs written to: {args.out_dir}")
    print("")
    print("Primary time-restricted result: vortioxetine vs escitalopram")
    vte = pairwise[pairwise["comparator"].eq("escitalopram")]
    if not vte.empty:
        r = vte.iloc[0]
        print(f"  {int(r['a_index_with_SD'])}/{int(r['n_index'])} vs {int(r['c_comparator_with_SD'])}/{int(r['n_comparator'])}; ROR {r['ror']:.3f} ({r['ror025']:.3f}-{r['ror975']:.3f}); {r['direction']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
