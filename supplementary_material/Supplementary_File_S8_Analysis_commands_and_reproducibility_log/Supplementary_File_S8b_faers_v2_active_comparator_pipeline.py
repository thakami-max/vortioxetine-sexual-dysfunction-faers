#!/usr/bin/env python3
"""
Focused V2 FAERS/AEMS active-comparator pipeline for vortioxetine-associated sexual dysfunction.

Designed for FDA FAERS/AEMS quarterly ASCII ZIP files.

Primary use case:
  - index drug: vortioxetine
  - active comparators: escitalopram, paroxetine, sertraline, duloxetine, venlafaxine, bupropion, mirtazapine
  - indication restriction: depression-related indications using INDI linked by PRIMARYID + DRUG_SEQ
  - exposure definition: primary suspect (PS) reports
  - monotherapy-like restriction: exclude reports containing >1 selected antidepressant coded PS or SS
  - stratification: overall, sex-stratified, optional country/report-source strata

Outputs are case/non-case pairwise Reporting Odds Ratios (RORs). These are signal-detection
estimates, not incidence, prevalence, relative risk, or proof of causality.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None


# ----------------------------- basic utilities ----------------------------- #

def log(msg: str) -> None:
    print(msg, flush=True)


def normalize_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    s = str(value).upper()
    s = s.replace("Æ", "AE").replace("Œ", "OE")
    s = re.sub(r"[^A-Z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def normalize_series(s: pd.Series) -> pd.Series:
    return (
        s.astype("string")
        .fillna("")
        .str.upper()
        .str.replace("Æ", "AE", regex=False)
        .str.replace("Œ", "OE", regex=False)
        .str.replace(r"[^A-Z0-9]+", " ", regex=True)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )


def detect_sep_from_zip_member(zip_path: Path, member: str) -> str:
    with zipfile.ZipFile(zip_path) as z:
        with z.open(member) as fh:
            sample = fh.read(4096).decode("latin1", errors="ignore")
    counts = {"$": sample.count("$"), "\t": sample.count("\t"), ",": sample.count(","), "|": sample.count("|")}
    return max(counts, key=counts.get) if max(counts.values()) > 0 else "$"


def iter_zip_members(raw_dir: Path, table_prefix: str) -> Iterable[Tuple[Path, str]]:
    table_prefix = table_prefix.upper()
    for zip_path in sorted(raw_dir.glob("*.zip")):
        try:
            with zipfile.ZipFile(zip_path) as z:
                for name in z.namelist():
                    base = Path(name).name.upper()
                    if base.startswith(table_prefix) and base.endswith((".TXT", ".CSV", ".ASC")):
                        yield zip_path, name
        except zipfile.BadZipFile:
            log(f"WARNING: skipping corrupted/non-ZIP file: {zip_path}")


def iter_deleted_members(raw_dir: Path) -> Iterable[Tuple[Path, str]]:
    for zip_path in sorted(raw_dir.glob("*.zip")):
        try:
            with zipfile.ZipFile(zip_path) as z:
                for name in z.namelist():
                    base = Path(name).name.upper()
                    if (base.startswith("DELETE") or base.startswith("DELETED")) and base.endswith((".TXT", ".CSV", ".ASC")):
                        yield zip_path, name
        except zipfile.BadZipFile:
            continue


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [str(c).strip().lower() for c in df.columns]
    return df


def read_member(zip_path: Path, member: str, chunksize: Optional[int] = None) -> Iterable[pd.DataFrame]:
    sep = detect_sep_from_zip_member(zip_path, member)
    with zipfile.ZipFile(zip_path) as z:
        with z.open(member) as fh:
            reader = pd.read_csv(
                fh,
                sep=sep,
                dtype=str,
                encoding="latin1",
                low_memory=False,
                on_bad_lines="skip",
                chunksize=chunksize,
            )
            if chunksize:
                for chunk in reader:
                    yield standardize_columns(chunk)
            else:
                yield standardize_columns(reader)


def find_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    cols = set(df.columns)
    for c in candidates:
        c = c.lower()
        if c in cols:
            return c
    return None


def require_col(df: pd.DataFrame, candidates: Sequence[str], table_name: str) -> str:
    col = find_col(df, candidates)
    if col is None:
        raise ValueError(f"Could not find any of columns {candidates} in {table_name}; found {list(df.columns)[:25]}")
    return col


def parse_faers_date(s: pd.Series) -> pd.Series:
    raw = s.astype("string").str.replace(r"[^0-9]", "", regex=True)
    dt = pd.to_datetime(raw, format="%Y%m%d", errors="coerce")
    fallback = pd.to_datetime(s, errors="coerce")
    return dt.fillna(fallback)


def convert_age_to_years(age: pd.Series, age_cod: pd.Series) -> pd.Series:
    val = pd.to_numeric(age, errors="coerce")
    unit = age_cod.astype("string").str.upper().str.strip()
    years = pd.Series(np.nan, index=age.index, dtype="float")
    years = years.mask(unit.isin(["YR", "Y", "YEAR", "YEARS"]), val)
    years = years.mask(unit.isin(["MON", "MO", "MONTH", "MONTHS"]), val / 12)
    years = years.mask(unit.isin(["WK", "W", "WEEK", "WEEKS"]), val / 52.1429)
    years = years.mask(unit.isin(["DY", "D", "DAY", "DAYS"]), val / 365.25)
    years = years.mask(unit.isin(["HR", "H", "HOUR", "HOURS"]), val / 8766)
    years = years.mask(unit.isin(["DEC", "DECADE", "DECADES"]), val * 10)
    return years


def age_group_from_years(years: pd.Series) -> pd.Series:
    out = pd.Series("Unknown", index=years.index, dtype="string")
    out = out.mask((years >= 0) & (years < 18), "0-17")
    out = out.mask((years >= 18) & (years < 45), "18-44")
    out = out.mask((years >= 45) & (years < 65), "45-64")
    out = out.mask(years >= 65, "65+")
    return out


def standardize_sex(sex: pd.Series) -> pd.Series:
    raw = sex.astype("string").str.upper().str.strip()
    out = pd.Series("Unknown", index=sex.index, dtype="string")
    out = out.mask(raw.isin(["M", "MALE"]), "Male")
    out = out.mask(raw.isin(["F", "FEMALE"]), "Female")
    return out


def standardize_reporter_type(occp: pd.Series) -> pd.Series:
    raw = occp.astype("string").str.upper().str.strip()
    out = pd.Series("Unknown/Other", index=occp.index, dtype="string")
    # FAERS occp_cod examples: MD, PH, RN, HP, OT, CN
    out = out.mask(raw.isin(["MD", "PH", "RN", "HP", "DENTIST", "LAWYER"]), "Healthcare professional")
    out = out.mask(raw.isin(["CN", "CONSUMER", "PATIENT"]), "Consumer/non-HCP")
    return out


def standardize_country(country: pd.Series) -> pd.Series:
    raw = country.astype("string").fillna("").str.upper().str.strip()
    raw = raw.str.replace(r"[^A-Z]+", " ", regex=True).str.strip()
    out = raw.mask(raw.isin(["UNITED STATES", "UNITED STATES OF AMERICA", "USA", "U S A", "US"]), "US")
    out = out.mask(out == "", "Unknown")
    return out


# ----------------------------- data loading ----------------------------- #

def load_deleted_identifiers(raw_dir: Path, chunksize: int) -> Tuple[Set[str], Set[str]]:
    primaryids: Set[str] = set()
    caseids: Set[str] = set()
    members = list(iter_deleted_members(raw_dir))
    if not members:
        log("No DELETE/DELETED files detected in ZIP members.")
        return primaryids, caseids

    log(f"Reading deleted-case files from {len(members)} ZIP member(s)...")
    for zip_path, member in members:
        log(f"  DELETE: {zip_path.name}::{Path(member).name}")
        try:
            for df in read_member(zip_path, member, chunksize=chunksize):
                if df.empty:
                    continue
                primary_col = find_col(df, ["primaryid", "isr"])
                case_col = find_col(df, ["caseid", "case", "case_num"])
                if primary_col:
                    primaryids.update(df[primary_col].astype("string").str.strip().dropna().tolist())
                if case_col:
                    caseids.update(df[case_col].astype("string").str.strip().dropna().tolist())
        except Exception as e:
            log(f"WARNING: could not parse deleted file {zip_path.name}::{member}: {e}")
    log(f"Deleted identifiers collected: {len(primaryids):,} primaryids; {len(caseids):,} caseids")
    return primaryids, caseids


def load_demo(raw_dir: Path, chunksize: int, remove_deleted: bool = True) -> pd.DataFrame:
    deleted_primaryids: Set[str] = set()
    deleted_caseids: Set[str] = set()
    if remove_deleted:
        deleted_primaryids, deleted_caseids = load_deleted_identifiers(raw_dir, chunksize=chunksize)

    frames: List[pd.DataFrame] = []
    members = list(iter_zip_members(raw_dir, "DEMO"))
    if not members:
        raise FileNotFoundError(f"No DEMO*.TXT/CSV files found inside ZIPs in {raw_dir}")
    log(f"Reading DEMO tables from {len(members)} ZIP member(s)...")
    for zip_path, member in members:
        log(f"  DEMO: {zip_path.name}::{Path(member).name}")
        for df in read_member(zip_path, member, chunksize=chunksize):
            if df.empty:
                continue
            primary_col = require_col(df, ["primaryid", "isr"], "DEMO")
            case_col = require_col(df, ["caseid", "case", "case_num"], "DEMO")
            fda_dt_col = find_col(df, ["fda_dt", "receiptdate", "mfr_dt"])
            optional = [
                fda_dt_col,
                find_col(df, ["age"]),
                find_col(df, ["age_cod"]),
                find_col(df, ["sex"]),
                find_col(df, ["occp_cod"]),
                find_col(df, ["reporter_country", "reporter_cntry"]),
                find_col(df, ["occr_country", "occur_country"]),
                find_col(df, ["rept_cod"]),
                find_col(df, ["rept_dt"]),
                find_col(df, ["mfr_dt"]),
                find_col(df, ["init_fda_dt"]),
            ]
            keep_cols = [primary_col, case_col] + [c for c in optional if c and c not in [primary_col, case_col]]
            sub = df[keep_cols].copy()
            rename = {primary_col: "primaryid", case_col: "caseid"}
            if fda_dt_col:
                rename[fda_dt_col] = "fda_dt"
            sub = sub.rename(columns=rename)
            frames.append(sub)

    demo = pd.concat(frames, ignore_index=True)
    demo["primaryid"] = demo["primaryid"].astype("string").str.strip()
    demo["caseid"] = demo["caseid"].astype("string").str.strip()
    demo = demo[(demo["primaryid"].notna()) & (demo["primaryid"] != "") & (demo["caseid"].notna()) & (demo["caseid"] != "")].copy()

    if remove_deleted and (deleted_primaryids or deleted_caseids):
        before = len(demo)
        demo = demo[~demo["primaryid"].isin(deleted_primaryids) & ~demo["caseid"].isin(deleted_caseids)].copy()
        log(f"Removed deleted reports from DEMO rows: {before - len(demo):,}")

    if "fda_dt" not in demo.columns:
        demo["fda_dt"] = pd.NA
    demo["fda_dt_parsed"] = parse_faers_date(demo["fda_dt"])
    demo["primaryid_num"] = pd.to_numeric(demo["primaryid"].str.extract(r"(\d+)")[0], errors="coerce")

    log(f"DEMO rows before case-level deduplication: {len(demo):,}")
    demo = demo.sort_values(["caseid", "fda_dt_parsed", "primaryid_num"], na_position="first")
    deduped = demo.groupby("caseid", as_index=False, sort=False).tail(1).copy()
    log(f"DEMO rows after case-level deduplication: {len(deduped):,}")

    if "age" in deduped.columns and "age_cod" in deduped.columns:
        deduped["age_years"] = convert_age_to_years(deduped["age"], deduped["age_cod"])
        deduped["age_group"] = age_group_from_years(deduped["age_years"])
    else:
        deduped["age_years"] = np.nan
        deduped["age_group"] = "Unknown"

    if "sex" in deduped.columns:
        deduped["sex_std"] = standardize_sex(deduped["sex"])
    else:
        deduped["sex_std"] = "Unknown"

    if "occp_cod" in deduped.columns:
        deduped["reporter_type_std"] = standardize_reporter_type(deduped["occp_cod"])
    else:
        deduped["reporter_type_std"] = "Unknown/Other"

    country_col = find_col(deduped, ["reporter_country", "reporter_cntry", "occr_country", "occur_country"])
    if country_col:
        deduped["country_std"] = standardize_country(deduped[country_col])
    else:
        deduped["country_std"] = "Unknown"

    deduped["is_us"] = deduped["country_std"].eq("US")
    return deduped.reset_index(drop=True)


# ----------------------------- exposure matching ----------------------------- #

@dataclass(frozen=True)
class DrugPattern:
    drug_class: str
    generic: str
    synonym: str
    analysis_role: str
    pattern: re.Pattern


def load_drug_patterns(path: Path, include_sensitivity: bool = False) -> Tuple[List[DrugPattern], pd.DataFrame]:
    cfg = pd.read_csv(path, dtype=str).fillna("")
    cfg.columns = [c.strip().lower() for c in cfg.columns]
    required = {"drug_class", "generic", "synonyms", "include_primary"}
    missing = required - set(cfg.columns)
    if missing:
        raise ValueError(f"Drug dictionary missing columns: {missing}")
    if not include_sensitivity:
        cfg = cfg[cfg["include_primary"].astype(str).str.strip().isin(["1", "true", "TRUE", "yes", "YES"])].copy()
    if "analysis_role" not in cfg.columns:
        cfg["analysis_role"] = "comparator"
    patterns: List[DrugPattern] = []
    seen: Set[Tuple[str, str, str]] = set()
    for _, row in cfg.iterrows():
        generic = normalize_text(row["generic"]).lower()
        drug_class = str(row["drug_class"]).strip()
        role = str(row.get("analysis_role", "comparator")).strip()
        syns = [row["generic"]] + str(row["synonyms"]).split(";")
        for syn in syns:
            norm = normalize_text(syn)
            if not norm:
                continue
            key = (drug_class, generic, norm)
            if key in seen:
                continue
            seen.add(key)
            pat = re.compile(r"(?<![A-Z0-9])" + re.escape(norm).replace(r"\ ", r"\s+") + r"(?![A-Z0-9])")
            patterns.append(DrugPattern(drug_class=drug_class, generic=generic, synonym=norm, analysis_role=role, pattern=pat))
    cfg["generic"] = cfg["generic"].map(lambda x: normalize_text(x).lower())
    return patterns, cfg


def extract_selected_drug_exposures(
    raw_dir: Path,
    valid_primaryids: Set[str],
    drug_patterns: List[DrugPattern],
    scan_role_codes: Sequence[str],
    chunksize: int,
) -> pd.DataFrame:
    members = list(iter_zip_members(raw_dir, "DRUG"))
    if not members:
        raise FileNotFoundError(f"No DRUG*.TXT/CSV files found inside ZIPs in {raw_dir}")

    role_set = {r.strip().upper() for r in scan_role_codes if r.strip()}
    records: List[pd.DataFrame] = []
    log(f"Scanning DRUG tables from {len(members)} ZIP member(s) for selected antidepressant exposures...")
    for zip_path, member in members:
        log(f"  DRUG: {zip_path.name}::{Path(member).name}")
        for df in read_member(zip_path, member, chunksize=chunksize):
            if df.empty:
                continue
            primary_col = require_col(df, ["primaryid", "isr"], "DRUG")
            seq_col = find_col(df, ["drug_seq", "drugseq", "seq_num", "drug_seq_num"])
            drug_col = find_col(df, ["drugname", "drug_name"])
            prod_ai_col = find_col(df, ["prod_ai", "active_ing", "active_ingredient"])
            role_col = find_col(df, ["role_cod", "role"])
            if not drug_col and not prod_ai_col:
                continue

            work = df.copy()
            work["primaryid"] = work[primary_col].astype("string").str.strip()
            work = work[work["primaryid"].isin(valid_primaryids)].copy()
            if work.empty:
                continue

            if role_col:
                work["role_cod"] = work[role_col].astype("string").str.upper().str.strip()
                if role_set:
                    work = work[work["role_cod"].isin(role_set)].copy()
                    if work.empty:
                        continue
            else:
                work["role_cod"] = ""

            if seq_col:
                work["drug_seq"] = work[seq_col].astype("string").str.strip()
            else:
                work["drug_seq"] = ""

            text_parts = []
            if drug_col:
                text_parts.append(normalize_series(work[drug_col]))
            if prod_ai_col:
                text_parts.append(normalize_series(work[prod_ai_col]))
            drug_text = text_parts[0] if text_parts else pd.Series("", index=work.index)
            for t in text_parts[1:]:
                drug_text = (drug_text + " " + t).str.strip()

            base = pd.DataFrame({
                "primaryid": work["primaryid"],
                "drug_seq": work["drug_seq"],
                "role_cod": work["role_cod"],
                "drug_text": drug_text,
            })
            matched_frames: List[pd.DataFrame] = []
            for patt in drug_patterns:
                mask = base["drug_text"].str.contains(patt.pattern, regex=True, na=False)
                if mask.any():
                    m = base.loc[mask, ["primaryid", "drug_seq", "role_cod", "drug_text"]].copy()
                    m["generic"] = patt.generic
                    m["drug_class"] = patt.drug_class
                    m["analysis_role"] = patt.analysis_role
                    m["matched_synonym"] = patt.synonym
                    matched_frames.append(m)
            if matched_frames:
                records.append(pd.concat(matched_frames, ignore_index=True))
    if not records:
        return pd.DataFrame(columns=["primaryid", "drug_seq", "role_cod", "generic", "drug_class", "analysis_role", "matched_synonym", "drug_text"])
    exposures = pd.concat(records, ignore_index=True)
    exposures = exposures.drop_duplicates(subset=["primaryid", "drug_seq", "role_cod", "generic", "drug_class"])
    log(f"Matched selected-antidepressant exposure rows: {len(exposures):,}; unique exposed cases: {exposures['primaryid'].nunique():,}")
    return exposures


# ----------------------------- indication and event matching ----------------------------- #

def load_indication_terms(path: Path, include_sensitivity: bool = False) -> pd.DataFrame:
    cfg = pd.read_csv(path, dtype=str).fillna("")
    cfg.columns = [c.strip().lower() for c in cfg.columns]
    required = {"indication", "include_primary"}
    missing = required - set(cfg.columns)
    if missing:
        raise ValueError(f"Indication dictionary missing columns: {missing}")
    if not include_sensitivity:
        cfg = cfg[cfg["include_primary"].astype(str).str.strip().isin(["1", "true", "TRUE", "yes", "YES"])].copy()
    cfg["indication_norm"] = cfg["indication"].map(normalize_text)
    cfg = cfg[cfg["indication_norm"] != ""].drop_duplicates("indication_norm")
    return cfg


def extract_indications(
    raw_dir: Path,
    valid_primaryids: Set[str],
    indication_cfg: pd.DataFrame,
    chunksize: int,
) -> pd.DataFrame:
    members = list(iter_zip_members(raw_dir, "INDI"))
    if not members:
        log("WARNING: No INDI files found; indication-restricted analysis cannot be performed.")
        return pd.DataFrame(columns=["primaryid", "drug_seq", "indi_pt", "indi_norm", "depression_indication"])

    terms = set(indication_cfg["indication_norm"].tolist())
    term_regexes = [(term, re.compile(r"(?<![A-Z0-9])" + re.escape(term).replace(r"\ ", r"\s+") + r"(?![A-Z0-9])")) for term in terms]
    records: List[pd.DataFrame] = []

    log(f"Scanning INDI tables from {len(members)} ZIP member(s) for depression-related indications...")
    for zip_path, member in members:
        log(f"  INDI: {zip_path.name}::{Path(member).name}")
        for df in read_member(zip_path, member, chunksize=chunksize):
            if df.empty:
                continue
            primary_col = require_col(df, ["primaryid", "isr"], "INDI")
            seq_col = find_col(df, ["indi_drug_seq", "drug_seq", "drugseq", "seq_num"])
            pt_col = require_col(df, ["indi_pt", "pt", "indication", "indi"], "INDI")
            work = df.copy()
            work["primaryid"] = work[primary_col].astype("string").str.strip()
            work = work[work["primaryid"].isin(valid_primaryids)].copy()
            if work.empty:
                continue
            work["drug_seq"] = work[seq_col].astype("string").str.strip() if seq_col else ""
            work["indi_pt"] = work[pt_col].astype("string").fillna("")
            work["indi_norm"] = normalize_series(work[pt_col])
            # Exact match or token-boundary phrase match
            exact = work["indi_norm"].isin(terms)
            phrase = pd.Series(False, index=work.index)
            if not exact.all():
                for _, pat in term_regexes:
                    phrase = phrase | work["indi_norm"].str.contains(pat, regex=True, na=False)
            work["depression_indication"] = exact | phrase
            keep = work[work["depression_indication"]].copy()
            if not keep.empty:
                records.append(keep[["primaryid", "drug_seq", "indi_pt", "indi_norm", "depression_indication"]].drop_duplicates())

    if not records:
        log("No depression-related indications matched in INDI.")
        return pd.DataFrame(columns=["primaryid", "drug_seq", "indi_pt", "indi_norm", "depression_indication"])
    out = pd.concat(records, ignore_index=True).drop_duplicates()
    log(f"Matched depression-related indication rows: {len(out):,}; unique cases: {out['primaryid'].nunique():,}")
    return out


def load_event_terms(path: Path, include_sensitivity: bool = False) -> pd.DataFrame:
    cfg = pd.read_csv(path, dtype=str).fillna("")
    cfg.columns = [c.strip().lower() for c in cfg.columns]
    required = {"domain", "pt", "include_primary"}
    missing = required - set(cfg.columns)
    if missing:
        raise ValueError(f"Event dictionary missing columns: {missing}")
    if not include_sensitivity:
        cfg = cfg[cfg["include_primary"].astype(str).str.strip().isin(["1", "true", "TRUE", "yes", "YES"])].copy()
    cfg["pt_norm"] = cfg["pt"].map(normalize_text)
    cfg = cfg[cfg["pt_norm"] != ""].drop_duplicates("pt_norm")
    return cfg


def extract_sd_events(raw_dir: Path, valid_primaryids: Set[str], term_cfg: pd.DataFrame, chunksize: int) -> pd.DataFrame:
    members = list(iter_zip_members(raw_dir, "REAC"))
    if not members:
        raise FileNotFoundError(f"No REAC*.TXT/CSV files found inside ZIPs in {raw_dir}")
    term_map = term_cfg.set_index("pt_norm").to_dict("index")
    term_set = set(term_map.keys())
    records: List[pd.DataFrame] = []
    log(f"Scanning REAC tables from {len(members)} ZIP member(s) for sexual dysfunction terms...")
    for zip_path, member in members:
        log(f"  REAC: {zip_path.name}::{Path(member).name}")
        for df in read_member(zip_path, member, chunksize=chunksize):
            if df.empty:
                continue
            primary_col = require_col(df, ["primaryid", "isr"], "REAC")
            pt_col = require_col(df, ["pt", "event_pt", "reac_pt"], "REAC")
            work = df.copy()
            work["primaryid"] = work[primary_col].astype("string").str.strip()
            work = work[work["primaryid"].isin(valid_primaryids)].copy()
            if work.empty:
                continue
            work["pt_norm"] = normalize_series(work[pt_col])
            m = work[work["pt_norm"].isin(term_set)].copy()
            if m.empty:
                continue
            m["pt_reported"] = m[pt_col]
            m["pt_dictionary"] = m["pt_norm"].map(lambda x: term_map[x]["pt"])
            m["domain"] = m["pt_norm"].map(lambda x: term_map[x]["domain"])
            records.append(m[["primaryid", "pt_reported", "pt_dictionary", "pt_norm", "domain"]].drop_duplicates())
    if not records:
        return pd.DataFrame(columns=["primaryid", "pt_reported", "pt_dictionary", "pt_norm", "domain"])
    events = pd.concat(records, ignore_index=True).drop_duplicates()
    log(f"Matched sexual dysfunction event rows: {len(events):,}; unique event cases: {events['primaryid'].nunique():,}")
    return events


# ----------------------------- cohort building ----------------------------- #

def apply_indication_restriction(
    exposures: pd.DataFrame,
    depression_indications: pd.DataFrame,
    mode: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if mode == "none":
        out = exposures.copy()
        out["depression_indication"] = pd.NA
        out["matched_indication_terms"] = ""
        summary = pd.DataFrame([{"restriction": "none", "rows_in": len(exposures), "rows_out": len(out), "excluded_rows": 0}])
        return out, summary

    if depression_indications.empty:
        log("WARNING: indication restriction requested but no depression indications were matched.")
        out = exposures.iloc[0:0].copy()
        summary = pd.DataFrame([{"restriction": mode, "rows_in": len(exposures), "rows_out": 0, "excluded_rows": len(exposures)}])
        return out, summary

    indi = depression_indications.copy()
    agg = indi.groupby(["primaryid", "drug_seq"]).agg(
        depression_indication=("depression_indication", "max"),
        matched_indication_terms=("indi_pt", lambda x: "; ".join(sorted(set(map(str, x))))),
    ).reset_index()
    merged = exposures.merge(agg, on=["primaryid", "drug_seq"], how="left")
    merged["depression_indication"] = merged["depression_indication"].fillna(False).astype(bool)
    merged["matched_indication_terms"] = merged["matched_indication_terms"].fillna("")
    out = merged[merged["depression_indication"]].copy()
    log(f"Indication restriction '{mode}': kept {len(out):,}/{len(exposures):,} exposure rows ({out['primaryid'].nunique():,} unique cases).")
    summary = pd.DataFrame([{
        "restriction": mode,
        "rows_in": len(exposures),
        "rows_out": len(out),
        "excluded_rows": len(exposures) - len(out),
        "unique_cases_out": out["primaryid"].nunique() if not out.empty else 0,
    }])
    return out, summary


def apply_monotherapy_like(
    exposures_primary: pd.DataFrame,
    exposures_for_exclusion: pd.DataFrame,
    monotherapy_like: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not monotherapy_like:
        summary = pd.DataFrame([{
            "monotherapy_like": False,
            "rows_in": len(exposures_primary),
            "rows_out": len(exposures_primary),
            "excluded_rows": 0,
            "unique_cases_out": exposures_primary["primaryid"].nunique() if not exposures_primary.empty else 0,
        }])
        return exposures_primary.copy(), summary

    if exposures_for_exclusion.empty or exposures_primary.empty:
        return exposures_primary.copy(), pd.DataFrame()

    gen_by_case = exposures_for_exclusion.groupby("primaryid")["generic"].apply(lambda x: sorted(set(map(str, x)))).reset_index()
    gen_by_case["n_selected_antidepressants_psss"] = gen_by_case["generic"].map(len)
    gen_by_case["only_selected_antidepressant"] = gen_by_case["generic"].map(lambda xs: xs[0] if len(xs) == 1 else "")
    gen_by_case = gen_by_case.drop(columns=["generic"])
    merged = exposures_primary.merge(gen_by_case, on="primaryid", how="left")
    before = len(merged)
    keep = merged[(merged["n_selected_antidepressants_psss"] == 1) & (merged["only_selected_antidepressant"] == merged["generic"])].copy()
    log(f"Monotherapy-like restriction: kept {len(keep):,}/{before:,} exposure rows ({keep['primaryid'].nunique():,} unique cases).")
    summary = pd.DataFrame([{
        "monotherapy_like": True,
        "rows_in": before,
        "rows_out": len(keep),
        "excluded_rows": before - len(keep),
        "unique_cases_out": keep["primaryid"].nunique() if not keep.empty else 0,
    }])
    return keep, summary


def ensure_single_final_exposure(exposures: pd.DataFrame) -> Tuple[pd.DataFrame, int]:
    if exposures.empty:
        return exposures.copy(), 0
    n_by_case = exposures.groupby("primaryid")["generic"].nunique()
    multi_cases = set(n_by_case[n_by_case > 1].index.astype(str))
    out = exposures[~exposures["primaryid"].astype(str).isin(multi_cases)].copy()
    return out, len(multi_cases)


# ----------------------------- statistics ----------------------------- #

def contingency_counts(index_cases: Set[str], comparator_cases: Set[str], outcome_cases: Set[str], universe: Optional[Set[str]] = None) -> Tuple[int, int, int, int]:
    if universe is not None:
        index_cases = index_cases & universe
        comparator_cases = comparator_cases & universe
        outcome_cases = outcome_cases & universe
    # Remove overlap if present to maintain a clean pairwise contrast.
    overlap = index_cases & comparator_cases
    if overlap:
        index_cases = index_cases - overlap
        comparator_cases = comparator_cases - overlap
    a = len(index_cases & outcome_cases)
    b = len(index_cases - outcome_cases)
    c = len(comparator_cases & outcome_cases)
    d = len(comparator_cases - outcome_cases)
    return a, b, c, d


def disproportionality_stats(a: int, b: int, c: int, d: int) -> Dict[str, float]:
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
    return {"ror": ror, "ror025": lo, "ror975": hi, "prr": prr, "chi2": chi2, "haldane_correction": corrected}


def case_sets_by_generic(exposures: pd.DataFrame) -> Dict[str, Set[str]]:
    return {g: set(x["primaryid"].astype(str)) for g, x in exposures.groupby("generic")}


def case_sets_by_class(exposures: pd.DataFrame, exclude_index: str) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = {}
    for cls, x in exposures[exposures["generic"] != exclude_index].groupby("drug_class"):
        out[f"pooled_{normalize_text(cls).lower().replace(' ', '_')}"] = set(x["primaryid"].astype(str))
    return out


def pairwise_table(
    exposures: pd.DataFrame,
    outcome_cases: Set[str],
    demo: pd.DataFrame,
    index_drug: str,
    comparator_names: Sequence[str],
    min_cases: int,
    stratum_name: str = "Overall",
    stratum_value: str = "Overall",
    stratum_cases: Optional[Set[str]] = None,
    comparator_case_sets: Optional[Dict[str, Set[str]]] = None,
) -> pd.DataFrame:
    index_drug = index_drug.lower()
    generic_sets = case_sets_by_generic(exposures)
    if comparator_case_sets is None:
        comparator_case_sets = {name: generic_sets.get(name.lower(), set()) for name in comparator_names}
    index_cases = generic_sets.get(index_drug, set())
    rows: List[Dict[str, object]] = []

    for comp_label, comp_cases in comparator_case_sets.items():
        a, b, c, d = contingency_counts(index_cases, comp_cases, outcome_cases, universe=stratum_cases)
        stats = disproportionality_stats(a, b, c, d)
        rows.append({
            "stratum": stratum_name,
            "stratum_value": stratum_value,
            "index_drug": index_drug,
            "comparator": comp_label,
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
        })
    return pd.DataFrame(rows)


def pairwise_by_outcome_set(
    exposures: pd.DataFrame,
    outcome_sets: Dict[str, Set[str]],
    index_drug: str,
    comparator_names: Sequence[str],
    min_cases: int,
    outcome_label_col: str,
) -> pd.DataFrame:
    generic_sets = case_sets_by_generic(exposures)
    index_cases = generic_sets.get(index_drug.lower(), set())
    rows: List[Dict[str, object]] = []
    for outcome_label, outcome_cases in outcome_sets.items():
        for comp in comparator_names:
            comp_cases = generic_sets.get(comp.lower(), set())
            a, b, c, d = contingency_counts(index_cases, comp_cases, outcome_cases, universe=None)
            if max(a, c) < min_cases:
                continue
            stats = disproportionality_stats(a, b, c, d)
            rows.append({
                outcome_label_col: outcome_label,
                "index_drug": index_drug.lower(),
                "comparator": comp.lower(),
                "a_index_with_event": a,
                "b_index_without_event": b,
                "c_comparator_with_event": c,
                "d_comparator_without_event": d,
                "n_index": a + b,
                "n_comparator": c + d,
                **stats,
                "signal_higher_vortioxetine_lower95_gt1": bool(a >= min_cases and stats["ror025"] > 1),
                "signal_lower_vortioxetine_upper95_lt1": bool(c >= min_cases and stats["ror975"] < 1),
            })
    return pd.DataFrame(rows)


# ----------------------------- outputs and plots ----------------------------- #

def make_case_level_flags(demo: pd.DataFrame, exposures: pd.DataFrame, events: pd.DataFrame) -> pd.DataFrame:
    case_df = demo.copy()
    case_df["primaryid"] = case_df["primaryid"].astype(str)
    outcome_cases = set(events["primaryid"].astype(str)) if not events.empty else set()
    exposed_cases = set(exposures["primaryid"].astype(str)) if not exposures.empty else set()
    case_df["sd_any"] = case_df["primaryid"].isin(outcome_cases)
    case_df["selected_antidepressant_exposure"] = case_df["primaryid"].isin(exposed_cases)
    if not exposures.empty:
        generics = exposures.groupby("primaryid")["generic"].apply(lambda x: ";".join(sorted(set(map(str, x))))).rename("selected_generics")
        classes = exposures.groupby("primaryid")["drug_class"].apply(lambda x: ";".join(sorted(set(map(str, x))))).rename("selected_classes")
        case_df = case_df.merge(generics, left_on="primaryid", right_index=True, how="left")
        case_df = case_df.merge(classes, left_on="primaryid", right_index=True, how="left")
    if not events.empty:
        domains = events.groupby("primaryid")["domain"].apply(lambda x: ";".join(sorted(set(map(str, x))))).rename("sd_domains")
        pts = events.groupby("primaryid")["pt_dictionary"].apply(lambda x: ";".join(sorted(set(map(str, x))))).rename("sd_terms")
        case_df = case_df.merge(domains, left_on="primaryid", right_index=True, how="left")
        case_df = case_df.merge(pts, left_on="primaryid", right_index=True, how="left")
    return case_df


def save_pairwise_forest(results: pd.DataFrame, out_path: Path, title: str) -> None:
    if plt is None or results.empty:
        return
    df = results[(results["stratum"] == "Overall") & (results["stratum_value"] == "Overall")].copy()
    if df.empty:
        return
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["ror", "ror025", "ror975"])
    if df.empty:
        return
    df = df.sort_values("ror")
    y = np.arange(len(df))
    plt.figure(figsize=(9, max(4.5, 0.45 * len(df) + 2)))
    plt.errorbar(df["ror"], y, xerr=[df["ror"] - df["ror025"], df["ror975"] - df["ror"]], fmt="o", capsize=3)
    plt.axvline(1, linestyle="--")
    labels = df["comparator"].astype(str) + " (a=" + df["a_index_with_SD"].astype(str) + ", c=" + df["c_comparator_with_SD"].astype(str) + ")"
    plt.yticks(y, labels)
    plt.xscale("log")
    plt.xlabel("Pairwise ROR: vortioxetine vs comparator (log scale)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()


def write_summary(
    out_dir: Path,
    demo: pd.DataFrame,
    exposures_all: pd.DataFrame,
    exposures_final: pd.DataFrame,
    events: pd.DataFrame,
    args: argparse.Namespace,
    indication_summary: pd.DataFrame,
    monotherapy_summary: pd.DataFrame,
    n_multi_final_excluded: int,
) -> None:
    lines: List[str] = []
    lines.append("# Focused V2 FAERS/AEMS active-comparator pipeline summary")
    lines.append("")
    lines.append(f"Raw data directory: `{args.raw_dir}`")
    lines.append(f"Index drug: `{args.index_drug}`")
    lines.append(f"Comparators: `{', '.join(args.comparators)}`")
    lines.append(f"Role codes for exposure cohort: `{', '.join(args.role_codes)}`")
    lines.append(f"Role codes for monotherapy-like exclusion: `{', '.join(args.exclude_other_selected_roles)}`")
    lines.append(f"Indication mode: `{args.indication_mode}`")
    lines.append(f"Monotherapy-like restriction: `{args.monotherapy_like}`")
    lines.append(f"Remove deleted cases: `{not args.keep_deleted}`")
    lines.append("")
    lines.append("## Counts")
    lines.append(f"- Deduplicated FAERS/AEMS cases: {len(demo):,}")
    lines.append(f"- Selected-antidepressant exposure rows before restrictions: {len(exposures_all):,}")
    lines.append(f"- Selected-antidepressant exposure cases before restrictions: {exposures_all['primaryid'].nunique() if not exposures_all.empty else 0:,}")
    lines.append(f"- Final exposure rows after indication/role/monotherapy restrictions: {len(exposures_final):,}")
    lines.append(f"- Final exposure cases after restrictions: {exposures_final['primaryid'].nunique() if not exposures_final.empty else 0:,}")
    lines.append(f"- Cases excluded for >1 final selected exposure: {n_multi_final_excluded:,}")
    lines.append(f"- Sexual dysfunction event rows matched: {len(events):,}")
    lines.append(f"- Unique sexual dysfunction event cases: {events['primaryid'].nunique() if not events.empty else 0:,}")
    lines.append("")
    lines.append("## Interpretation note")
    lines.append("FAERS/AEMS is a spontaneous reporting system. Pairwise disproportionality results are signal-detection outputs, not incidence, prevalence, relative risk, or proof of causality.")
    lines.append("For pairwise RORs, values >1 indicate higher sexual-dysfunction reporting with vortioxetine relative to the comparator; values <1 indicate lower reporting with vortioxetine relative to the comparator.")
    lines.append("")
    lines.append("## Restriction summaries")
    if not indication_summary.empty:
        lines.append(indication_summary.to_string(index=False))
    if not monotherapy_summary.empty:
        lines.append("")
        lines.append(monotherapy_summary.to_string(index=False))
    (out_dir / "pipeline_v2_summary.md").write_text("\n".join(lines), encoding="utf-8")


# ----------------------------- main ----------------------------- #

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Focused V2 FAERS/AEMS active-comparator analysis for vortioxetine and sexual dysfunction")
    p.add_argument("--raw-dir", type=Path, default=Path("data_raw"))
    p.add_argument("--out-dir", type=Path, default=Path("outputs_v2_active_comparator"))
    p.add_argument("--drug-config", type=Path, default=Path("config/antidepressants_selected_comparators.csv"))
    p.add_argument("--event-config", type=Path, default=Path("config/sexual_dysfunction_terms_narrow_primary.csv"))
    p.add_argument("--indication-config", type=Path, default=Path("config/indications_depression.csv"))
    p.add_argument("--index-drug", default="vortioxetine")
    p.add_argument("--comparators", nargs="+", default=["escitalopram", "paroxetine", "sertraline", "duloxetine", "venlafaxine", "bupropion", "mirtazapine"])
    p.add_argument("--role-codes", nargs="+", default=["PS"], help="Exposure role codes for primary cohort; recommended primary: PS")
    p.add_argument("--exclude-other-selected-roles", nargs="+", default=["PS", "SS"], help="Selected-drug roles used to define monotherapy-like exclusion")
    p.add_argument("--indication-mode", choices=["depression", "none"], default="depression")
    p.add_argument("--monotherapy-like", action="store_true", help="Exclude cases with >1 selected antidepressant coded PS/SS")
    p.add_argument("--include-sensitivity-indications", action="store_true")
    p.add_argument("--include-sensitivity-terms", action="store_true")
    p.add_argument("--keep-deleted", action="store_true", help="Do not attempt to remove DELETE/DELETED case files")
    p.add_argument("--chunksize", type=int, default=250000)
    p.add_argument("--min-cases", type=int, default=3)
    p.add_argument("--run-sex-strata", action="store_true")
    p.add_argument("--run-source-strata", action="store_true", help="Add US-only and HCP-only/consumer strata")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    args.raw_dir = args.raw_dir.resolve()
    args.out_dir = args.out_dir.resolve()
    for sub in ["processed", "results", "figures"]:
        (args.out_dir / sub).mkdir(parents=True, exist_ok=True)

    log("=== Focused V2 FAERS/AEMS active-comparator analysis ===")
    log(f"Index drug: {args.index_drug}")
    log(f"Comparators: {', '.join(args.comparators)}")

    demo = load_demo(args.raw_dir, chunksize=args.chunksize, remove_deleted=not args.keep_deleted)
    valid_primaryids = set(demo["primaryid"].astype(str))
    demo.to_csv(args.out_dir / "processed" / "demo_deduplicated.csv", index=False)

    drug_patterns, drug_cfg = load_drug_patterns(args.drug_config, include_sensitivity=False)
    log(f"Loaded {len(drug_patterns)} selected drug-name patterns.")
    scan_roles = sorted(set([r.upper() for r in args.role_codes] + [r.upper() for r in args.exclude_other_selected_roles]))
    exposures_all = extract_selected_drug_exposures(args.raw_dir, valid_primaryids, drug_patterns, scan_roles, chunksize=args.chunksize)
    exposures_all.to_csv(args.out_dir / "processed" / "selected_antidepressant_exposures_all_scanned_roles.csv", index=False)

    # Build primary-role exposure rows first.
    role_set = {r.upper() for r in args.role_codes}
    exposures_primary_role = exposures_all[exposures_all["role_cod"].str.upper().isin(role_set)].copy()
    exposures_primary_role.to_csv(args.out_dir / "processed" / "selected_antidepressant_exposures_primary_roles_pre_indication.csv", index=False)

    indication_summary = pd.DataFrame()
    if args.indication_mode == "depression":
        indication_cfg = load_indication_terms(args.indication_config, include_sensitivity=args.include_sensitivity_indications)
        log(f"Loaded {len(indication_cfg)} depression-related indication terms.")
        depression_indications = extract_indications(args.raw_dir, valid_primaryids, indication_cfg, chunksize=args.chunksize)
        depression_indications.to_csv(args.out_dir / "processed" / "depression_indication_matches.csv", index=False)
        exposures_indication, indication_summary = apply_indication_restriction(exposures_primary_role, depression_indications, args.indication_mode)
    else:
        depression_indications = pd.DataFrame()
        exposures_indication, indication_summary = apply_indication_restriction(exposures_primary_role, depression_indications, "none")

    # Monotherapy-like exclusion uses all selected antidepressants in PS/SS roles, independent of indication.
    mono_role_set = {r.upper() for r in args.exclude_other_selected_roles}
    exposures_for_mono = exposures_all[exposures_all["role_cod"].str.upper().isin(mono_role_set)].copy()
    exposures_mono, monotherapy_summary = apply_monotherapy_like(exposures_indication, exposures_for_mono, monotherapy_like=args.monotherapy_like)

    # De-duplicate case-generic rows and remove any residual multiple final selected exposures.
    exposures_final = exposures_mono.drop_duplicates(subset=["primaryid", "generic", "drug_class"]).copy()
    exposures_final, n_multi_final_excluded = ensure_single_final_exposure(exposures_final)
    exposures_final.to_csv(args.out_dir / "processed" / "selected_antidepressant_exposures_final.csv", index=False)

    term_cfg = load_event_terms(args.event_config, include_sensitivity=args.include_sensitivity_terms)
    log(f"Loaded {len(term_cfg)} sexual dysfunction event terms.")
    events = extract_sd_events(args.raw_dir, valid_primaryids, term_cfg, chunksize=args.chunksize)
    events.to_csv(args.out_dir / "processed" / "sexual_dysfunction_event_matches.csv", index=False)

    case_flags = make_case_level_flags(demo, exposures_final, events)
    case_flags.to_csv(args.out_dir / "processed" / "case_level_flags_v2.csv", index=False)

    outcome_cases = set(events["primaryid"].astype(str)) if not events.empty else set()

    # Overall pairwise drug-vs-drug comparisons.
    overall = pairwise_table(
        exposures_final,
        outcome_cases,
        demo,
        index_drug=args.index_drug,
        comparator_names=args.comparators,
        min_cases=args.min_cases,
        stratum_name="Overall",
        stratum_value="Overall",
        stratum_cases=None,
    )
    overall.to_csv(args.out_dir / "results" / "pairwise_drug_comparisons_overall.csv", index=False)

    # Pooled SSRI/SNRI context.
    class_sets = case_sets_by_class(exposures_final, exclude_index=args.index_drug.lower())
    pooled_keep = {k: v for k, v in class_sets.items() if k in ["pooled_ssri", "pooled_snri"]}
    pooled = pairwise_table(
        exposures_final,
        outcome_cases,
        demo,
        index_drug=args.index_drug,
        comparator_names=list(pooled_keep.keys()),
        min_cases=args.min_cases,
        comparator_case_sets=pooled_keep,
    )
    pooled.to_csv(args.out_dir / "results" / "pairwise_pooled_class_comparisons_overall.csv", index=False)

    # Stratified analyses.
    strata_tables = [overall.assign(stratum="Overall", stratum_value="Overall")]
    if args.run_sex_strata:
        for sex in ["Male", "Female", "Unknown"]:
            cases = set(demo.loc[demo["sex_std"].astype(str).eq(sex), "primaryid"].astype(str))
            st = pairwise_table(exposures_final, outcome_cases, demo, args.index_drug, args.comparators, args.min_cases, "sex_std", sex, cases)
            strata_tables.append(st)

    if args.run_source_strata:
        source_specs = {
            "US_only": set(demo.loc[demo["is_us"].fillna(False), "primaryid"].astype(str)),
            "non_US_or_unknown": set(demo.loc[~demo["is_us"].fillna(False), "primaryid"].astype(str)),
            "HCP_only": set(demo.loc[demo["reporter_type_std"].eq("Healthcare professional"), "primaryid"].astype(str)),
            "consumer_non_HCP": set(demo.loc[demo["reporter_type_std"].eq("Consumer/non-HCP"), "primaryid"].astype(str)),
        }
        for label, cases in source_specs.items():
            st = pairwise_table(exposures_final, outcome_cases, demo, args.index_drug, args.comparators, args.min_cases, "source_country", label, cases)
            strata_tables.append(st)

    stratified = pd.concat(strata_tables, ignore_index=True) if strata_tables else overall
    stratified.to_csv(args.out_dir / "results" / "pairwise_drug_comparisons_stratified.csv", index=False)

    # Domain-specific and PT-specific pairwise results.
    if not events.empty:
        domain_sets = {domain: set(g["primaryid"].astype(str)) for domain, g in events.groupby("domain")}
        domain_results = pairwise_by_outcome_set(exposures_final, domain_sets, args.index_drug, args.comparators, args.min_cases, "domain")
        domain_results.to_csv(args.out_dir / "results" / "pairwise_domain_specific_comparisons.csv", index=False)

        pt_sets = {pt: set(g["primaryid"].astype(str)) for pt, g in events.groupby("pt_dictionary")}
        pt_results = pairwise_by_outcome_set(exposures_final, pt_sets, args.index_drug, args.comparators, args.min_cases, "pt")
        pt_results.to_csv(args.out_dir / "results" / "pairwise_pt_specific_comparisons.csv", index=False)

    # Counts
    if not exposures_final.empty:
        exp_counts = exposures_final.groupby(["drug_class", "generic"]).agg(
            unique_cases=("primaryid", "nunique"),
            matched_rows=("primaryid", "size"),
        ).reset_index().sort_values(["drug_class", "generic"])
        exp_counts.to_csv(args.out_dir / "results" / "selected_exposure_counts_final.csv", index=False)

        sd_by_generic = exposures_final.assign(sd_any=exposures_final["primaryid"].astype(str).isin(outcome_cases)).groupby(["drug_class", "generic"]).agg(
            n_exposed=("primaryid", "nunique"),
            sd_cases=("sd_any", "sum"),
        ).reset_index()
        sd_by_generic["sd_reporting_proportion"] = sd_by_generic["sd_cases"] / sd_by_generic["n_exposed"]
        sd_by_generic.to_csv(args.out_dir / "results" / "selected_exposure_sd_counts_final.csv", index=False)

    if not events.empty:
        event_counts = events.groupby(["domain", "pt_dictionary"]).agg(
            unique_cases=("primaryid", "nunique"),
            matched_rows=("primaryid", "size"),
        ).reset_index().sort_values(["unique_cases", "pt_dictionary"], ascending=[False, True])
        event_counts.to_csv(args.out_dir / "results" / "event_counts.csv", index=False)

    if not args.no_plots:
        save_pairwise_forest(overall, args.out_dir / "figures" / "forest_pairwise_vortioxetine_vs_comparators.png", "Vortioxetine vs selected antidepressant comparators")

    write_summary(args.out_dir, demo, exposures_all, exposures_final, events, args, indication_summary, monotherapy_summary, n_multi_final_excluded)
    log(f"Done. Outputs written to: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
