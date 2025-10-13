#!/usr/bin/env python3
import pandas as pd
import argparse

PALINDROMIC = {("A","T"), ("T","A"), ("C","G"), ("G","C")}
VALID_ALLELES = {"A","C","G","T"}

def prepare_for_ldsc(path_in, path_out):
    df = pd.read_csv(path_in, sep="\t")

    rename_map = {
        "Effect": "A1",
        "Non_Effect": "A2",
        "Beta": "BETA",
        "PVAL": "P",
        "P": "P",
        "Effect_allele_freq": "MAF",
        "IMPINFO": "INFO",
        "ID": "SNP"
    }
    df = df.rename(columns={k:v for k,v in rename_map.items() if k in df.columns})
    required = ["SNP", "A1", "A2", "BETA", "SE", "P"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["A1"] = df["A1"].astype(str).str.upper()
    df["A2"] = df["A2"].astype(str).str.upper()
    df = df[df["A1"].isin(VALID_ALLELES) & df["A2"].isin(VALID_ALLELES)].copy()

    pal = df.apply(lambda r: (r["A1"], r["A2"]) in PALINDROMIC, axis=1)
    if "MAF" in df.columns:
        df["MAF"] = pd.to_numeric(df["MAF"], errors="coerce")
        ambiguous = pal & df["MAF"].between(0.4, 0.6, inclusive="both")
        df = df[~ambiguous].copy()
    else:
        df = df[~pal].copy()

    if "INFO" in df.columns:
        df["INFO"] = pd.to_numeric(df["INFO"], errors="coerce")
        df = df[df["INFO"] >= 0.9]

    if "MAF" in df.columns:
        df = df[df["MAF"] >= 0.01]

    keep_cols = required + [c for c in ["MAF", "INFO"] if c in df.columns]
    df = df[keep_cols].dropna(subset=required)
    print(f"Number of SNPs after QC: {len(df):,}")

    df.to_csv(path_out, sep="\t", index=False, compression="infer")

def n_to_neff():
    print(f"\n AD Neff = {round(4 / (1/21982 + 1/41944), 0):,}")
    print(f" SZ Neff = {round(4 / (1/74776 + 1/101023), 0):,}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", required=True)
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()
    n_to_neff()
    prepare_for_ldsc(args.infile, args.outfile)
