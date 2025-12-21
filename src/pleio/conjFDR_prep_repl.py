#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def norm_cols(df):
    df.columns = [str(c).strip() for c in df.columns]
    return df

def read_sumstats(path):
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df = norm_cols(df)
    if {"SNP","CHR","POS"}.issubset(df.columns):
        out = pd.DataFrame()
        out["SNP"] = df["SNP"]
        out["CHR"] = pd.to_numeric(df["CHR"], errors="coerce").astype("Int64")
        out["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
        for c in ["A1","A2","P","BETA","SE","N"]:
            if c in df.columns:
                out[c] = df[c]
        out["build"] = "GRCh37"
        return out

    if {"rs","chr","hg19","ref","alt","allele1","beta","se","pvalue"}.issubset(df.columns):
        out = pd.DataFrame()
        out["SNP"] = df["rs"].astype(str)
        out["CHR"] = pd.to_numeric(df["chr"], errors="coerce").astype("Int64")
        out["POS"] = pd.to_numeric(df["hg19"], errors="coerce").astype("Int64")
        a1 = df["allele1"].astype(str).str.upper()
        ref = df["ref"].astype(str).str.upper()
        alt = df["alt"].astype(str).str.upper()
        out["A1"] = a1
        out["A2"] = ref.where(a1 != ref, alt)
        out["BETA"] = pd.to_numeric(df["beta"], errors="coerce")
        out["SE"] = pd.to_numeric(df["se"], errors="coerce")
        out["P"] = pd.to_numeric(df["pvalue"], errors="coerce")
        out["build"] = "GRCh37"
        return out
    raise SystemExit(f"Unrecognized columns in {path}: {list(df.columns)[:25]}")

def remove_mhc(df):
    m = ~((df["CHR"] == 6) & (df["POS"] >= 25119106) & (df["POS"] <= 33854733))
    return df[m].copy()

def prepare_trait(df, prefix):
    need = ["SNP","A1","A2","P","CHR","POS"]
    extra = [c for c in ["BETA","SE","N"] if c in df.columns]
    df = df[need + extra].copy()
    df["A1"] = df["A1"].astype(str).str.upper()
    df["A2"] = df["A2"].astype(str).str.upper()
    df = remove_mhc(df)
    df = df.drop_duplicates(subset=["SNP"])
    df = df[df["CHR"].between(1, 22)]
    df = df[df["POS"].notna()]
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df[df["P"].between(0, 1)]
    df = df.rename(columns={c: f"{c}_{prefix}" for c in df.columns if c not in ["SNP","CHR","POS"]})
    return df

def harmonise(t1, t2, p1, p2):
    m = t1.merge(t2, on=["SNP","CHR","POS"], how="inner")
    same = (m[f"A1_{p1}"] == m[f"A1_{p2}"]) & (m[f"A2_{p1}"] == m[f"A2_{p2}"])
    flip = (m[f"A1_{p1}"] == m[f"A2_{p2}"]) & (m[f"A2_{p1}"] == m[f"A1_{p2}"])
    m = m[same | flip].copy()

    b2 = f"BETA_{p2}"
    if b2 in m.columns:
        m.loc[flip, b2] = -pd.to_numeric(m.loc[flip, b2], errors="coerce")

    m.loc[flip, f"A1_{p2}"] = m.loc[flip, f"A1_{p1}"]
    m.loc[flip, f"A2_{p2}"] = m.loc[flip, f"A2_{p1}"]
    m["build"] = "GRCh37"
    return m

def build_cfdr_input(m, p1, p2):
    out = m[["SNP", f"P_{p1}", f"P_{p2}"]].copy()
    out = out.rename(columns={f"P_{p1}":"p1", f"P_{p2}":"p2"})
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--trait1", required=True)
    p.add_argument("--trait2", required=True)
    p.add_argument("--prefix1", required=True)
    p.add_argument("--prefix2", required=True)
    p.add_argument("--out_prefix", required=True)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    t1 = read_sumstats(Path(args.trait1))
    t2 = read_sumstats(Path(args.trait2))
    t1p = prepare_trait(t1, args.prefix1)
    t2p = prepare_trait(t2, args.prefix2)
    m = harmonise(t1p, t2p, args.prefix1, args.prefix2)
    cfdr = build_cfdr_input(m, args.prefix1, args.prefix2)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    m.to_csv(out_dir / f"{args.out_prefix}.harmonised_{args.out_prefix}.tsv", sep="\t", index=False)
    cfdr.to_csv(out_dir / f"{args.out_prefix}.cfdr_input_{args.out_prefix}.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
