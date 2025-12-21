#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

def load_tsv(p):
    return pd.read_csv(p, sep="\t", dtype=str, low_memory=False)

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

def exclude_mhc(d, cc="CHR", pc="BP"):
    c = d[cc]
    p = d[pc]
    mhc = (c == 6) & (p >= 25_000_000) & (p <= 34_000_000)
    return d[~mhc].copy()

def add_freqs(d):
    if "Effect_allele_freq" in d.columns:
        eaf = pd.to_numeric(d["Effect_allele_freq"], errors="coerce")
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

def filter_maf(d):
    return d[d["MAF"] >= 0.01].copy()

def rename_for_ldsc(d):
    m = {
        "CHR": "CHR",
        "BP": "POS",
        "SNP": "SNP",
        "Effect": "A1",
        "Non_Effect": "A2",
        "Beta": "BETA",
        "SE": "SE",
        "P": "P",
    }
    return d.rename(columns=m)

def dropna_required(d):
    req = ["SNP","A1","A2","FRQ","BETA","SE","P","CHR","POS"]
    return d.dropna(subset=[c for c in req if c in d.columns]).copy()

def write_ldsc(d, p):
    cols = ["SNP","A1","A2","FRQ","N","BETA","SE","P","CHR","POS"]
    if "N" not in d.columns:
        d["N"] = np.nan
    p.parent.mkdir(parents=True, exist_ok=True)
    d[cols].to_csv(p, sep="\t", index=False)

def count(lab, d):
    print(f"{lab}: {len(d):,}")

def parse_args():
    ap = argparse.ArgumentParser(description="QC + LDSC format for AD (Kunkle 2019)")
    ap.add_argument("--in",  dest="inp",  default="Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.tsv")
    ap.add_argument("--out", dest="out", default="Data/AD/post-qc/Kunkle_etal_2019_IGAP_Summary_statistics_published_ldsc_ready.tsv")
    return ap.parse_args()

if __name__ == "__main__":
    a = parse_args()
    src = Path(a.inp)
    out = Path(a.out)
    d = load_tsv(src)
    count("loaded", d)
    d = to_num(d, ["CHR","BP","Beta","SE","P","Effect_allele_freq","FRQ"])
    d = exclude_mhc(d, "CHR", "BP")
    count("after_exclude_MHC", d)
    d = keep_snps(d, "Effect", "Non_Effect")
    count("after_remove_indels", d)
    d = drop_palindromes(d, "Effect", "Non_Effect")
    count("after_remove_palindromes", d)
    d = add_freqs(d)
    d = filter_maf(d)
    count("after_maf", d)
    d = rename_for_ldsc(d)
    d = dropna_required(d)
    count("after_dropna_required", d)
    write_ldsc(d, out)
    print("wrote:", out)
