#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path("/Users/c24102394/Desktop/PhD/AD_SCZ_AGE")
CONJFDR_DIR = BASE / "outputs" / "conjFDR" / "hits"
OUT_DIR = CONJFDR_DIR / "shared"
THR = 0.05

def print_n(df: pd.DataFrame, label: str):
    print(f"Number of hits for {label}: {len(df):,}")

def print_n_shared(df: pd.DataFrame, df2: pd.DataFrame, label: str):
    counter = 0
    for _, row1 in df.iterrows():
        for _, row2 in df2.iterrows():
            if row1["SNP"] == row2["SNP"]:
                counter += 1
    print(f"Number of shared hits between {label}: {counter:,}")

def merge_shared_snps(df1: pd.DataFrame, df2: pd.DataFrame):
    return df1.merge(df2, on="SNP", how="inner")

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ad_long_path = CONJFDR_DIR / "AD_LONG_cfdr_results_mapped.tsv"
    scz_long_path = CONJFDR_DIR / "SCZ_LONG_cfdr_results_mapped.tsv"
    ad_long_all = pd.read_csv(ad_long_path, sep="\t", dtype={"SNP": str})
    scz_long_all = pd.read_csv(scz_long_path, sep="\t", dtype={"SNP": str})
    ad_long_hits = ad_long_all[ad_long_all["conj_fdr"] < THR].copy()
    scz_long_hits = scz_long_all[scz_long_all["conj_fdr"] < THR].copy()
    print_n(ad_long_hits, "AD-LONG conjFDR<0.05")
    print_n(scz_long_hits, "SCZ-LONG conjFDR<0.05")
    print_n_shared(ad_long_hits, scz_long_hits, "AD-LONG and SCZ-LONG (conjFDR<0.05)")
    shared = merge_shared_snps(ad_long_hits, scz_long_hits)
    print(f"Merged shared variants: {len(shared):,}")
    out_file = OUT_DIR / "AD_LONG_SCZ_LONG_shared_conjFDR0.05.tsv"
    shared.to_csv(out_file, sep="\t", index=False)
    print(f"Saved merged shared hits to: {out_file}")

if __name__ == "__main__":
    main()