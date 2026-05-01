#!/usr/bin/env python3
"""
Concomitant drug/condition sensitivity analysis for the
focused V2 FAERS/AEMS active-comparator study.

This script is designed to run inside the existing faers_antidepressant_sd_pipeline
folder, after the focused V2 primary analysis has been completed.

It:
  * loads final V2 selected antidepressant exposure cases;
  * loads sexual-dysfunction case flags from the V2 processed output;
  * scans raw FAERS/AEMS DRUG tables for concomitant drug groups associated with
    sexual dysfunction;
  * scans raw FAERS/AEMS INDI tables for urologic/endocrine and related condition
    terms;
  * computes pairwise vortioxetine-vs-comparator RORs after exclusions and in
    strata defined by confounder flags.

Important interpretation note:
  These sensitivity analyses reduce some observable confounding in public FAERS
  extracts, but they do not eliminate unmeasured confounding. PDE5-inhibitor
  co-reporting is flagged separately because it may indicate baseline erectile
  dysfunction or treatment of the reported outcome.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# ---------- utilities ----------

def norm_col(x: object) -> str:
    return str(x).strip().lower().replace(" ", "_").replace("-", "_")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [norm_col(c) for c in df.columns]
    return df


def find_col(df: pd.DataFrame, candidates: Sequence[str], required: bool = True) -> Optional[str]:
    cols = set(df.columns)
    for c in candidates:
        cc = norm_col(c)
        if cc in cols:
            return cc
    if required:
        raise KeyError(f"None of the candidate columns found: {candidates}; available columns include {list(df.columns)[:25]}")
    return None


def clean_text_series(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.lower().str.replace(r"\s+", " ", regex=True)


def compile_term_regex(terms: Sequence[str]) -> re.Pattern:
    # Word-ish boundary: avoid matching morphine inside apomorphine, etc.
    cleaned = [str(t).strip().lower() for t in terms if str(t).strip()]
    cleaned = sorted(set(cleaned), key=len, reverse=True)
    escaped = [re.escape(t) for t in cleaned]
    if not escaped:
        # Regex that never matches
        return re.compile(r"a^", flags=re.I)
    pattern = r"(?<![a-z0-9])(?:" + "|".join(escaped) + r")(?![a-z0-9])"
    return re.compile(pattern, flags=re.I)


def read_csv_flexible(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return normalize_columns(df)


def list_zip_members(raw_dir: Path, table_prefix: str) -> List[Tuple[Path, str]]:
    pairs: List[Tuple[Path, str]] = []
    for zpath in sorted(raw_dir.glob("*.zip")):
        try:
            with zipfile.ZipFile(zpath) as zf:
                for name in zf.namelist():
                    base = os.path.basename(name).lower()
                    if base.startswith(table_prefix.lower()) and base.endswith((".txt", ".csv")):
                        pairs.append((zpath, name))
        except zipfile.BadZipFile:
            print(f"WARNING: Skipping bad zip file: {zpath}", file=sys.stderr)
    return pairs


def iter_faers_table(raw_dir: Path, table_prefix: str, chunksize: int) -> Iterable[Tuple[Path, str, pd.DataFrame]]:
    members = list_zip_members(raw_dir, table_prefix)
    if not members:
        print(f"WARNING: no {table_prefix} tables found in {raw_dir}", file=sys.stderr)
    for zpath, member in members:
        print(f"  {table_prefix.upper()}: {zpath.name}::{os.path.basename(member)}", flush=True)
        with zipfile.ZipFile(zpath) as zf:
            with zf.open(member) as fh:
                try:
                    reader = pd.read_csv(
                        fh,
                        sep="$",
                        dtype=str,
                        chunksize=chunksize,
                        encoding="latin-1",
                        on_bad_lines="skip",
                        low_memory=False,
                    )
                    for chunk in reader:
                        yield zpath, member, normalize_columns(chunk)
                except Exception as e:
                    print(f"WARNING: failed with C engine for {zpath.name}::{member}: {e}; retrying python engine", file=sys.stderr)
                    fh.seek(0)
                    reader = pd.read_csv(
                        fh,
                        sep="$",
                        dtype=str,
                        chunksize=chunksize,
                        encoding="latin-1",
                        on_bad_lines="skip",
                        engine="python",
                    )
                    for chunk in reader:
                        yield zpath, member, normalize_columns(chunk)


# ---------- loading V2 outputs ----------

def load_final_exposures(v2_out_dir: Path) -> pd.DataFrame:
    candidates = [
        v2_out_dir / "processed" / "selected_antidepressant_exposures_final.csv",
        v2_out_dir / "processed" / "final_selected_exposures.csv",
        v2_out_dir / "processed" / "selected_exposures_final.csv",
    ]
    path = next((p for p in candidates if p.exists()), None)
    if path is None:
        raise FileNotFoundError(
            "Could not find final V2 exposure file. Expected one of:\n" + "\n".join(str(p) for p in candidates)
        )
    df = read_csv_flexible(path)
    pid = find_col(df, ["primaryid", "primary_id"])
    generic = find_col(df, ["generic", "matched_generic", "exposure", "drug", "index_or_comparator"], required=False)
    if generic is None:
        # Try to detect any column with generic/exposure in the name
        generic_candidates = [c for c in df.columns if "generic" in c or "exposure" in c]
        if generic_candidates:
            generic = generic_candidates[0]
        else:
            raise KeyError(f"Could not identify generic/exposure column in {path}; columns: {df.columns.tolist()}")
    out = df.copy()
    out["primaryid_std"] = out[pid].astype(str)
    out["generic_std"] = out[generic].astype(str).str.strip().str.lower()
    sex_col = find_col(out, ["sex", "sex_std", "gender"], required=False)
    if sex_col:
        out["sex_std"] = out[sex_col].fillna("Unknown").astype(str).str.strip()
    else:
        out["sex_std"] = "Unknown"
    # one row per primaryid/generic to avoid duplicate drug-name mappings
    out = out.drop_duplicates(["primaryid_std", "generic_std"]).reset_index(drop=True)
    return out


def load_sd_cases(v2_out_dir: Path) -> set:
    candidates = [
        v2_out_dir / "processed" / "sexual_dysfunction_event_matches.csv",
        v2_out_dir / "processed" / "sd_event_matches.csv",
        v2_out_dir / "processed" / "case_level_flags_v2.csv",
    ]
    # Prefer event matches, because column names are simpler.
    for path in candidates:
        if path.exists():
            df = read_csv_flexible(path)
            if "has_sd" in df.columns:
                pid = find_col(df, ["primaryid", "primary_id"])
                return set(df.loc[df["has_sd"].astype(str).str.lower().isin(["true", "1", "yes"]), pid].astype(str))
            pid = find_col(df, ["primaryid", "primary_id"], required=False)
            if pid:
                return set(df[pid].dropna().astype(str).unique())
    raise FileNotFoundError(
        "Could not locate sexual-dysfunction event or case flag file in V2 outputs."
    )


# ---------- scan confounders ----------

def load_term_config(path: Path) -> pd.DataFrame:
    df = read_csv_flexible(path)
    group_col = find_col(df, ["group"])
    term_col = find_col(df, ["term"])
    df = df[[group_col, term_col] + [c for c in df.columns if c not in {group_col, term_col}]].copy()
    df = df.rename(columns={group_col: "group", term_col: "term"})
    df["group"] = df["group"].astype(str).str.strip().str.lower()
    df["term"] = df["term"].astype(str).str.strip().str.lower()
    df = df[(df["group"] != "") & (df["term"] != "")]
    return df


def build_group_regexes(cfg: pd.DataFrame) -> Dict[str, re.Pattern]:
    return {g: compile_term_regex(sub["term"].tolist()) for g, sub in cfg.groupby("group")}


def scan_drug_confounders(raw_dir: Path, cfg: pd.DataFrame, chunksize: int, case_universe: Optional[set] = None) -> pd.DataFrame:
    regexes = build_group_regexes(cfg)
    rows: List[pd.DataFrame] = []
    print(f"Scanning DRUG tables for {len(regexes)} concomitant drug groups...")
    for _zpath, _member, chunk in iter_faers_table(raw_dir, "DRUG", chunksize):
        pid = find_col(chunk, ["primaryid", "primary_id"], required=False)
        if not pid:
            continue
        if case_universe is not None:
            chunk = chunk[chunk[pid].astype(str).isin(case_universe)]
            if chunk.empty:
                continue
        # Avoid sequence/role columns but use drug/product/active ingredient columns.
        text_cols = []
        for c in chunk.columns:
            cl = c.lower()
            if cl in {"primaryid", "caseid", "drug_seq", "role_cod", "val_vbm", "route", "dose_vbm", "cum_dose_chr", "cum_dose_unit", "dechal", "rechal"}:
                continue
            if any(k in cl for k in ["drug", "prod", "active", "ingredient", "medicinal"]):
                text_cols.append(c)
        if not text_cols:
            # fallback: scan all non-ID text columns
            text_cols = [c for c in chunk.columns if c not in {pid}]
        text = clean_text_series(chunk[text_cols].astype(str).agg(" ".join, axis=1))
        pid_series = chunk[pid].astype(str)
        for group, rgx in regexes.items():
            mask = text.str.contains(rgx)
            if mask.any():
                rows.append(pd.DataFrame({"primaryid_std": pid_series[mask].values, group: True}))
    if not rows:
        return pd.DataFrame(columns=["primaryid_std"] + list(regexes.keys()))
    flags = pd.concat(rows, ignore_index=True)
    # Collapse flags by primaryid.
    for group in regexes:
        if group not in flags.columns:
            flags[group] = False
    agg = flags.groupby("primaryid_std", as_index=False)[list(regexes.keys())].max()
    return agg


def scan_condition_confounders(raw_dir: Path, cfg: pd.DataFrame, chunksize: int, case_universe: Optional[set] = None) -> pd.DataFrame:
    regexes = build_group_regexes(cfg)
    rows: List[pd.DataFrame] = []
    print(f"Scanning INDI tables for {len(regexes)} condition groups...")
    for _zpath, _member, chunk in iter_faers_table(raw_dir, "INDI", chunksize):
        pid = find_col(chunk, ["primaryid", "primary_id"], required=False)
        if not pid:
            continue
        if case_universe is not None:
            chunk = chunk[chunk[pid].astype(str).isin(case_universe)]
            if chunk.empty:
                continue
        text_cols = []
        for c in chunk.columns:
            cl = c.lower()
            if cl in {"primaryid", "caseid", "drug_seq"}:
                continue
            if any(k in cl for k in ["indi", "pt", "indication", "term"]):
                text_cols.append(c)
        if not text_cols:
            text_cols = [c for c in chunk.columns if c not in {pid}]
        text = clean_text_series(chunk[text_cols].astype(str).agg(" ".join, axis=1))
        pid_series = chunk[pid].astype(str)
        for group, rgx in regexes.items():
            mask = text.str.contains(rgx)
            if mask.any():
                rows.append(pd.DataFrame({"primaryid_std": pid_series[mask].values, group: True}))
    if not rows:
        return pd.DataFrame(columns=["primaryid_std"] + list(regexes.keys()))
    flags = pd.concat(rows, ignore_index=True)
    for group in regexes:
        if group not in flags.columns:
            flags[group] = False
    agg = flags.groupby("primaryid_std", as_index=False)[list(regexes.keys())].max()
    return agg


# ---------- statistics ----------

def calc_stats(a: int, b: int, c: int, d: int) -> Dict[str, float]:
    # Haldane-Anscombe correction for zero cells only.
    aa, bb, cc, dd = map(float, [a, b, c, d])
    corrected = False
    if min(aa, bb, cc, dd) == 0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
        corrected = True
    ror = (aa / bb) / (cc / dd)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    lo = math.exp(math.log(ror) - 1.96 * se)
    hi = math.exp(math.log(ror) + 1.96 * se)
    # PRR and chi-square, mostly for robustness reporting.
    p1 = aa / (aa + bb)
    p0 = cc / (cc + dd)
    prr = p1 / p0 if p0 > 0 else np.nan
    n = aa + bb + cc + dd
    denom = (aa + cc) * (bb + dd) * (aa + bb) * (cc + dd)
    chi2 = n * (aa * dd - bb * cc) ** 2 / denom if denom > 0 else np.nan
    return {
        "ror": ror,
        "ror025": lo,
        "ror975": hi,
        "prr": prr,
        "chi2": chi2,
        "zero_cell_corrected": corrected,
    }


def pairwise_table(
    df: pd.DataFrame,
    index_drug: str,
    comparators: Sequence[str],
    min_cases: int,
    analysis_label: str,
) -> pd.DataFrame:
    rows: List[dict] = []
    index_drug = index_drug.lower()
    comps = [c.lower() for c in comparators]
    for comp in comps:
        sub = df[df["generic_std"].isin([index_drug, comp])].copy()
        idx = sub["generic_std"].eq(index_drug)
        cmp = sub["generic_std"].eq(comp)
        a = int((idx & sub["has_sd"]).sum())
        b = int((idx & ~sub["has_sd"]).sum())
        c = int((cmp & sub["has_sd"]).sum())
        d = int((cmp & ~sub["has_sd"]).sum())
        stats = calc_stats(a, b, c, d)
        row = {
            "analysis": analysis_label,
            "index_drug": index_drug,
            "comparator": comp,
            "a_index_with_SD": a,
            "b_index_without_SD": b,
            "n_index": a + b,
            "c_comparator_with_SD": c,
            "d_comparator_without_SD": d,
            "n_comparator": c + d,
            **stats,
        }
        row["signal_higher_vortioxetine_lower95_gt1"] = bool(a >= min_cases and stats["ror025"] > 1)
        row["signal_lower_vortioxetine_upper95_lt1"] = bool(a >= min_cases and stats["ror975"] < 1)
        row["prr_signal_higher_vortioxetine"] = bool(a >= min_cases and stats["prr"] >= 2 and stats["chi2"] >= 4)
        if stats["ror975"] < 1:
            row["direction"] = "lower_with_vortioxetine"
        elif stats["ror025"] > 1:
            row["direction"] = "higher_with_vortioxetine"
        else:
            row["direction"] = "no_clear_difference"
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_flags(df: pd.DataFrame, groups: Sequence[str]) -> pd.DataFrame:
    rows = []
    for drug, sub in df.groupby("generic_std"):
        base = {"generic": drug, "n_exposed": len(sub), "sd_cases": int(sub["has_sd"].sum())}
        for g in groups:
            base[f"{g}_cases"] = int(sub[g].sum()) if g in sub.columns else 0
            base[f"{g}_pct"] = (100 * base[f"{g}_cases"] / len(sub)) if len(sub) else 0
        rows.append(base)
    return pd.DataFrame(rows).sort_values("generic")


# ---------- main ----------

def main() -> int:
    p = argparse.ArgumentParser(description="Concomitant drug/condition sensitivity analysis for V2 FAERS/AEMS active-comparator study")
    p.add_argument("--raw-dir", default="data_raw")
    p.add_argument("--v2-out-dir", default="outputs_v2_primary_PS_depression_mono_narrow")
    p.add_argument("--out-dir", default="outputs_v2_confounder_sensitivity")
    p.add_argument("--drug-confounder-config", default="config/concomitant_drug_groups.csv")
    p.add_argument("--condition-confounder-config", default="config/concomitant_condition_terms.csv")
    p.add_argument("--index-drug", default="vortioxetine")
    p.add_argument("--comparators", nargs="+", default=["escitalopram", "paroxetine", "sertraline", "duloxetine", "venlafaxine", "bupropion", "mirtazapine"])
    p.add_argument("--chunksize", type=int, default=250000)
    p.add_argument("--min-cases", type=int, default=3)
    args = p.parse_args()

    raw_dir = Path(args.raw_dir)
    v2_out_dir = Path(args.v2_out_dir)
    out_dir = Path(args.out_dir)
    (out_dir / "processed").mkdir(parents=True, exist_ok=True)
    (out_dir / "results").mkdir(parents=True, exist_ok=True)

    print("=== FAERS/AEMS V2 concomitant drug/condition sensitivity analysis ===")
    print(f"Raw directory: {raw_dir}")
    print(f"V2 output directory: {v2_out_dir}")
    print(f"Output directory: {out_dir}")

    exposures = load_final_exposures(v2_out_dir)
    sd_cases = load_sd_cases(v2_out_dir)
    exposures["has_sd"] = exposures["primaryid_std"].isin(sd_cases)

    index_and_comps = [args.index_drug.lower()] + [c.lower() for c in args.comparators]
    exposures = exposures[exposures["generic_std"].isin(index_and_comps)].copy()
    case_universe = set(exposures["primaryid_std"].astype(str).unique())
    print(f"Loaded final V2 selected exposure cases: {len(exposures):,} rows; {len(case_universe):,} unique cases")
    print(f"Sexual-dysfunction cases in selected V2 cohort: {int(exposures['has_sd'].sum()):,}")

    drug_cfg = load_term_config(Path(args.drug_confounder_config))
    cond_cfg = load_term_config(Path(args.condition_confounder_config))

    drug_flags = scan_drug_confounders(raw_dir, drug_cfg, args.chunksize, case_universe=case_universe)
    cond_flags = scan_condition_confounders(raw_dir, cond_cfg, args.chunksize, case_universe=case_universe)

    drug_flags.to_csv(out_dir / "processed" / "concomitant_drug_flags_by_case.csv", index=False)
    cond_flags.to_csv(out_dir / "processed" / "condition_flags_by_case.csv", index=False)

    df = exposures.copy()
    for flags in [drug_flags, cond_flags]:
        if not flags.empty:
            df = df.merge(flags, on="primaryid_std", how="left")

    drug_groups = sorted(drug_cfg["group"].unique().tolist())
    cond_groups = sorted(cond_cfg["group"].unique().tolist())
    groups = drug_groups + cond_groups
    for g in groups:
        if g not in df.columns:
            df[g] = False
        df[g] = df[g].fillna(False).astype(bool)

    # Core clinically relevant groups. PDE5 is flagged separately because it may be an outcome-treatment proxy.
    pde5_group = "pde5_or_ed_treatment"
    core_groups = [g for g in groups if g != pde5_group]
    core_confounder_groups = [
        "antipsychotic",
        "five_ari_alpha_blocker",
        "opioid",
        "hormonal_therapy",
        "urologic_endocrine_condition",
    ]
    core_confounder_groups = [g for g in core_confounder_groups if g in groups]

    df["any_core_confounder"] = df[core_confounder_groups].any(axis=1) if core_confounder_groups else False
    df["any_confounder_excluding_pde5"] = df[core_groups].any(axis=1) if core_groups else False
    df["any_confounder_including_pde5"] = df[groups].any(axis=1) if groups else False

    df.to_csv(out_dir / "processed" / "final_exposures_with_confounder_flags.csv", index=False)

    # Summaries.
    flag_summary = summarize_flags(df, groups + ["any_core_confounder", "any_confounder_excluding_pde5", "any_confounder_including_pde5"])
    flag_summary.to_csv(out_dir / "results" / "confounder_flag_counts_by_drug.csv", index=False)

    all_results = []
    all_results.append(pairwise_table(df, args.index_drug, args.comparators, args.min_cases, "primary_v2_replicated"))
    all_results.append(pairwise_table(df[~df["any_core_confounder"]], args.index_drug, args.comparators, args.min_cases, "exclude_core_confounders_except_pde5"))
    all_results.append(pairwise_table(df[~df["any_confounder_excluding_pde5"]], args.index_drug, args.comparators, args.min_cases, "exclude_any_confounder_except_pde5"))
    all_results.append(pairwise_table(df[~df["any_confounder_including_pde5"]], args.index_drug, args.comparators, args.min_cases, "exclude_any_confounder_including_pde5"))

    # Exclude each group one at a time.
    for g in groups:
        all_results.append(pairwise_table(df[~df[g]], args.index_drug, args.comparators, args.min_cases, f"exclude_{g}"))

    # Stratify by yes/no for each core flag and individual group.
    strat_rows = []
    strat_flags = ["any_core_confounder", "any_confounder_excluding_pde5", "any_confounder_including_pde5"] + groups
    for flag in strat_flags:
        for val in [False, True]:
            sub = df[df[flag].eq(val)]
            if sub.empty:
                continue
            tab = pairwise_table(sub, args.index_drug, args.comparators, args.min_cases, f"stratified_{flag}_{'yes' if val else 'no'}")
            tab.insert(1, "stratum", flag)
            tab.insert(2, "stratum_value", "yes" if val else "no")
            strat_rows.append(tab)

    overall = pd.concat(all_results, ignore_index=True)
    stratified = pd.concat(strat_rows, ignore_index=True) if strat_rows else pd.DataFrame()
    overall.to_csv(out_dir / "results" / "pairwise_confounder_sensitivity_overall.csv", index=False)
    stratified.to_csv(out_dir / "results" / "pairwise_confounder_sensitivity_stratified.csv", index=False)

    # Compact manuscript-oriented table for primary comparator.
    primary_comp = "escitalopram"
    esc = overall[overall["comparator"].eq(primary_comp)].copy()
    esc.to_csv(out_dir / "results" / "primary_comparison_vortioxetine_vs_escitalopram_confounder_sensitivity.csv", index=False)

    # Summary markdown.
    summary = out_dir / "confounder_sensitivity_summary.md"
    with summary.open("w", encoding="utf-8") as f:
        f.write("# Concomitant drug/condition sensitivity analysis summary\n\n")
        f.write(f"V2 output directory: `{v2_out_dir}`\n\n")
        f.write(f"Final selected exposure rows loaded: {len(df):,}\n")
        f.write(f"Unique selected cases loaded: {df['primaryid_std'].nunique():,}\n")
        f.write(f"Sexual-dysfunction cases in selected cohort: {int(df['has_sd'].sum()):,}\n\n")
        f.write("## Confounder groups flagged\n\n")
        for g in groups:
            f.write(f"- {g}: {int(df[g].sum()):,} cases\n")
        f.write(f"- any core confounder excluding PDE5: {int(df['any_core_confounder'].sum()):,} cases\n")
        f.write(f"- any confounder excluding PDE5: {int(df['any_confounder_excluding_pde5'].sum()):,} cases\n")
        f.write(f"- any confounder including PDE5: {int(df['any_confounder_including_pde5'].sum()):,} cases\n\n")
        f.write("## Interpretation note\n\n")
        f.write("These analyses are sensitivity analyses based on observable drug and indication/reporting fields in public FAERS/AEMS extracts. They cannot remove unmeasured confounding. PDE5-inhibitor co-reporting was flagged separately because it may represent treatment of erectile dysfunction rather than baseline confounding.\n")

    print("Done. Outputs written to:", out_dir)
    print("Key files:")
    print(" -", out_dir / "confounder_sensitivity_summary.md")
    print(" -", out_dir / "results" / "pairwise_confounder_sensitivity_overall.csv")
    print(" -", out_dir / "results" / "pairwise_confounder_sensitivity_stratified.csv")
    print(" -", out_dir / "results" / "primary_comparison_vortioxetine_vs_escitalopram_confounder_sensitivity.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
