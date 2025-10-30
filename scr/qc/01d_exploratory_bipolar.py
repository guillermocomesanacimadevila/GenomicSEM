#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path
import numpy as np

def convert_to_tsv(file_path):
    df = pd.read_csv(file_path, sep=r"\s+", engine="python")
    df = df[["SNP","CHR","BP","A1","A2","P","SE","INFO","HRC_FRQ_A1","OR"]]
    df["BETA"] = np.log(df["OR"])
    return df

def to_num(d, cols):
    for c in cols:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d

def keep_snps(d, a1c, a2c):
    b = {"A","C","G","T"}
    a1 = d[a1c].astype(str).str.upper()
    a2 = d[a2c].astype(str).str.upper()
    snp = a1.isin(b) & a2.isin(b) & (a1.str.len()==1) & (a2.str.len()==1)
    gap = a1.str.contains("-", na=False) | a2.str.contains("-", na=False)
    return d[snp & ~gap].copy()

def drop_palindromes(d, a1c, a2c):
    a1 = d[a1c].astype(str).str.upper()
    a2 = d[a2c].astype(str).str.upper()
    pal = ((a1=="A")&(a2=="T")) | ((a1=="T")&(a2=="A")) | ((a1=="C")&(a2=="G")) | ((a1=="G")&(a2=="C"))
    return d[~pal].copy()

def exclude_regions(d, cc="CHR", pc="BP"):
    def inr(r,c,s,e):
        try:
            return int(r[cc]) == c and s <= int(r[pc]) <= e
        except:
            return False
    mhc  = d.apply(lambda r: inr(r, 6, 25_000_000, 34_000_000), axis=1)
    apoe = d.apply(lambda r: inr(r, 19, 45_116_911, 46_318_605), axis=1)
    return d[~(mhc | apoe)].copy()

def add_freqs(d):
    if "HRC_FRQ_A1" in d.columns:
        eaf = pd.to_numeric(d["HRC_FRQ_A1"], errors="coerce")
        d["FRQ"] = eaf
        d["MAF"] = np.minimum(eaf, 1 - eaf)
    elif "FRQ" in d.columns:
        maf = pd.to_numeric(d["FRQ"], errors="coerce")
        d["MAF"] = maf
        d["FRQ"] = maf
    else:
        d["FRQ"] = np.nan
        d["MAF"] = np.nan
    return d

def filter_maf_info(d):
    return d[(d["MAF"] >= 0.01) & (d["INFO"] >= 0.90)].copy()

def dropna_required(d):
    req = ["SNP","A1","A2","FRQ","BETA","SE","P","CHR","BP"]
    return d.dropna(subset=[c for c in req if c in d.columns]).copy()

def write_ldsc(d, p):
    cols = ["SNP","A1","A2","FRQ","N","BETA","SE","P","CHR","BP"]
    if "N" not in d.columns:
        d["N"] = np.nan
    p.parent.mkdir(parents=True, exist_ok=True)
    d[cols].to_csv(p, sep="\t", index=False)

def count(lab, d):
    print(f"{lab}: {len(d):,}")

def parse_args():
    ap = argparse.ArgumentParser(description="QC + LDSC format for BIP 2024 EUR no23andMe")
    ap.add_argument("--in",  dest="inp",  default="Data/BIP/bip2024_eur_no23andMe",
                    help="input BIP file (space-delimited)")
    ap.add_argument("--out", dest="out", default="Data/BIP/post-qc/bip2024_eur_no23andMe_ldsc_ready.tsv",
                    help="output LDSC-ready tsv")
    return ap.parse_args()

if __name__ == "__main__":
    a = parse_args()
    src = Path(a.inp)
    out = Path(a.out)
    d = convert_to_tsv(src)
    count("loaded", d)
    d = to_num(d, ["CHR","BP","BETA","SE","P","FRQ","INFO"])
    d = exclude_regions(d, "CHR", "BP")
    count("after_exclude_MHC_APOE", d)
    d = keep_snps(d, "A1", "A2")
    count("after_remove_indels", d)
    d = drop_palindromes(d, "A1", "A2")
    count("after_remove_palindromes", d)
    d = add_freqs(d)
    d = filter_maf_info(d)
    count("after_maf_info", d)
    d = dropna_required(d)
    count("after_dropna_required", d)
    write_ldsc(d, out)
    print("wrote:", out)
