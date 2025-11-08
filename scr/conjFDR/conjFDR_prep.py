#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def read_sumstats(path):
    df = pd.read_csv(path, sep="\t", dtype={"CHR": "Int64", "POS": "Int64"})
    df = df[df["CHR"].between(1, 22)]
    df["build"] = "GRCh37"
    return df

def remove_mhc(df):
    m = ~((df["CHR"] == 6) & (df["POS"] >= 25119106) & (df["POS"] <= 33854733))
    return df[m].copy()

def prepare_trait(df, prefix):
    cols = ["SNP", "A1", "A2", "P", "CHR", "POS"]
    extra = [c for c in ["BETA", "SE", "N"] if c in df.columns]
    df = df[cols + extra].copy()
    df = remove_mhc(df)
    df = df.drop_duplicates(subset=["SNP"])
    df = df[df["P"].between(0, 1)]
    df = df.rename(columns={c: f"{c}_{prefix}" for c in df.columns if c not in ["SNP", "CHR", "POS"]})
    return df

def harmonise(ad, scz):
    m = ad.merge(scz, on=["SNP", "CHR", "POS"], how="inner")
    same = (m["A1_AD"] == m["A1_SCZ"]) & (m["A2_AD"] == m["A2_SCZ"])
    flip = (m["A1_AD"] == m["A2_SCZ"]) & (m["A2_AD"] == m["A1_SCZ"])
    keep = same | flip
    m = m[keep].copy()
    if "BETA_SCZ" in m.columns:
        m.loc[flip, "BETA_SCZ"] = -m.loc[flip, "BETA_SCZ"]
    m.loc[flip, "A1_SCZ"] = m.loc[flip, "A1_AD"]
    m.loc[flip, "A2_SCZ"] = m.loc[flip, "A2_AD"]
    m["build"] = "GRCh37"
    return m

def build_cfdr_input(m):
    out = m[["SNP", "P_SCZ", "P_AD"]].copy()
    out = out.rename(columns={"P_SCZ": "p1", "P_AD": "p2"})
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ad", required=True)
    p.add_argument("--scz", required=True)
    p.add_argument("--out_prefix", required=True)
    args = p.parse_args()

    ad = read_sumstats(Path(args.ad))
    scz = read_sumstats(Path(args.scz))
    ad_prep = prepare_trait(ad, "AD")
    scz_prep = prepare_trait(scz, "SCZ")
    m = harmonise(ad_prep, scz_prep)
    cfdr = build_cfdr_input(m)

    out_dir = Path("/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/scr/conjFDR/data")
    out_dir.mkdir(parents=True, exist_ok=True)

    m.to_csv(out_dir / f"{args.out_prefix}.harmonised_AD_SCZ.tsv", sep="\t", index=False)
    cfdr.to_csv(out_dir / f"{args.out_prefix}.cfdr_input_AD_SCZ.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
