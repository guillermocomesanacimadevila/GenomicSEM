#!/usr/bin/env python3
import numpy as np
import pandas as pd
from pathlib import Path

def load(path):
    return pd.read_csv(path, sep="\t")

def rename_cols(df):
    df = df.drop(
        columns=[
            "snpid",
            "beta1_healthspan","se_healthspan",
            "beta1_lifespan","se_lifespan",
            "beta1_longevity","se_longevity",
        ],
        errors="ignore"
    )
    df = df.reset_index(drop=True).rename(
        columns={
            "rsid": "SNP",
            "chr": "CHR",
            "pos": "POS",
            "a1": "A1",
            "a0": "A2",
            "freq1": "EAF",
            "n": "N",
            "beta1": "BETA",
            "se": "SE",
            "info": "INFO",
            "p": "P",
        }
    )
    return df

def to_num(df):
    for c in ["CHR","POS","N","BETA","SE","P","EAF","INFO"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def exclude_regions(df):
    def inr(r,c,s,e):
        try: return int(r["CHR"])==c and s <= int(r["POS"]) <= e
        except: return False
    m_mhc  = df.apply(lambda r: inr(r,6,25_000_000,34_000_000), axis=1)
    m_apoe = df.apply(lambda r: inr(r,19,45_116_911,46_318_605), axis=1)
    return df[~(m_mhc | m_apoe)].copy()

def keep_snps(df):
    b = {"A","C","G","T"}
    a1 = df["A1"].astype(str).str.upper()
    a2 = df["A2"].astype(str).str.upper()
    snp = a1.isin(b) & a2.isin(b) & (a1.str.len()==1) & (a2.str.len()==1)
    gap = a1.str.contains("-", na=False) | a2.str.contains("-", na=False)
    return df[snp & ~gap].copy()

def drop_palindromes(df):
    a1 = df["A1"].astype(str).str.upper()
    a2 = df["A2"].astype(str).str.upper()
    pal = ((a1=="A")&(a2=="T"))|((a1=="T")&(a2=="A"))|((a1=="C")&(a2=="G"))|((a1=="G")&(a2=="C"))
    return df[~pal].copy()

def maf_info_filter(df):
    df["FRQ"] = df["EAF"]
    df["MAF"] = np.minimum(df["FRQ"], 1 - df["FRQ"])
    return df[(df["MAF"] >= 0.01) & (df["INFO"] >= 0.90)].copy()

def dropna_req(df):
    req = ["SNP","A1","A2","FRQ","N","BETA","SE","P","CHR","POS","INFO"]
    return df.dropna(subset=[c for c in req if c in df.columns]).copy()

def write_ldsc(df, path):
    cols = ["SNP","A1","A2","FRQ","N","BETA","SE","P","CHR","POS","INFO"]
    path.parent.mkdir(parents=True, exist_ok=True)
    df[cols].to_csv(path, sep="\t", index=False)

def count(label, df):
    print(f"{label}: {len(df):,}")

if __name__ == "__main__":
    src = Path("Data/AGE/timmers2020_healthspan_lifespan_longevity.tsv")
    out = Path("Data/AGE/post-qc/AGE_ldsc_ready.tsv")

    df = load(src)
    count("loaded", df)
    df = rename_cols(df)
    df = to_num(df)
    count("after_rename_numeric", df)
    df = exclude_regions(df)
    count("after_exclude_MHC_APOE", df)
    df = keep_snps(df)
    count("after_remove_indels", df)
    df = drop_palindromes(df)
    count("after_remove_palindromes", df)
    df = maf_info_filter(df)
    count("after_maf_info", df)
    df = dropna_req(df)
    count("after_dropna_required", df)
    write_ldsc(df, out)
    print("wrote:", out)
