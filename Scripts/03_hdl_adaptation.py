#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def load_sumstats(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)

def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)

def normalise_headers(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    rename_map = {
        "ID": "SNP",
        "Effect": "A1",
        "Non_Effect": "A2",
        "EA": "A1",
        "NEA": "A2",
        "Beta": "BETA",
        "PVAL": "P",
        "ImpINFO": "INFO",
        "IMPINFO": "INFO",
        "Effect_allele_freq": "MAF",
        "N_total": "N",
        "Ncases": "N_cases",
        "Ncontrols": "N_controls"
    }
    for k, v in rename_map.items():
        if k in d.columns and v not in d.columns:
            d = d.rename(columns={k: v})
    return d

def ensure_columns(d: pd.DataFrame) -> None:
    has_beta_se = ("BETA" in d.columns) and ("SE" in d.columns)
    has_z = ("Z" in d.columns)
    if not (has_beta_se or has_z):
        raise ValueError("Need either (BETA and SE) or Z column")

def uppercase_alleles(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    d["A1"] = d["A1"].astype(str).str.upper()
    d["A2"] = d["A2"].astype(str).str.upper()
    return d

def compute_z(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    if "BETA" in d.columns and "SE" in d.columns:
        b = pd.to_numeric(d["BETA"], errors="coerce")
        se = pd.to_numeric(d["SE"], errors="coerce")
        d["b"] = b
        d["se"] = se
        d["Z"] = b / se
    elif "Z" in d.columns:
        d["Z"] = pd.to_numeric(d["Z"], errors="coerce")
        if "b" not in d.columns:
            d["b"] = np.nan
        if "se" not in d.columns:
            d["se"] = np.nan
    return d

def add_constant_n(d: pd.DataFrame, n_value: float) -> pd.DataFrame:
    d = d.copy()
    d["N"] = float(n_value)
    return d

def select_hdl_columns(d: pd.DataFrame) -> pd.DataFrame:
    base_cols = ["SNP", "A1", "A2", "N", "Z", "b", "se"]
    for c in base_cols:
        if c not in d.columns:
            d[c] = np.nan
    extras = [c for c in ["MAF", "INFO"] if c in d.columns]
    return d[base_cols + extras]

def build_hdl_input(df: pd.DataFrame, n_value: float) -> pd.DataFrame:
    d = normalise_headers(df)
    ensure_columns(d)
    d = uppercase_alleles(d)
    d = compute_z(d)
    d = add_constant_n(d, n_value)
    d = select_hdl_columns(d)
    return d

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--n", type=float, required=True)
    return ap.parse_args()

def main() -> None:
    args = parse_args()
    df = load_sumstats(Path(args.infile))
    out = build_hdl_input(df, n_value=args.n)
    write_tsv(out, Path(args.outfile))
    print(f"Saved HDL-L file: {args.outfile} ({len(out):,} rows)")

if __name__ == "__main__":
    main()
