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

def exclude_mhc(d, cc="CHROM", pc="POS"):
    c = d[cc]
    p = d[pc]
    mhc = (c == 6) & (p >= 25_000_000) & (p <= 34_000_000)
    return d[~mhc].copy()

def compute_eaf(d):
    if {"FCAS","FCON"}.issubset(d.columns):
        fcas = pd.to_numeric(d["FCAS"], errors="coerce")
        fcon = pd.to_numeric(d["FCON"], errors="coerce")
        if {"NCAS","NCON"}.issubset(d.columns):
            ncas = pd.to_numeric(d["NCAS"], errors="coerce")
            ncon = pd.to_numeric(d["NCON"], errors="coerce")
            w = ncas + ncon
            eaf = np.where(w > 0, (fcas*ncas + fcon*ncon)/w, (fcas+fcon)/2.0)
        else:
            eaf = (fcas + fcon) / 2.0
        d["FRQ"] = eaf
        d["MAF"] = np.minimum(eaf, 1 - eaf)
    else:
        d["FRQ"] = np.nan
        d["MAF"] = np.nan
    return d

def add_info_and_n(d):
    if "IMPINFO" in d.columns and "INFO" not in d.columns:
        d = d.rename(columns={"IMPINFO": "INFO"})
    d["INFO"] = pd.to_numeric(d.get("INFO", np.nan), errors="coerce")
    if {"NCAS","NCON"}.issubset(d.columns):
        n = pd.to_numeric(d["NCAS"], errors="coerce") + pd.to_numeric(d["NCON"], errors="coerce")
        d["N"] = n
    elif "N" not in d.columns:
        d["N"] = np.nan
    return d

def filter_maf_info(d):
    return d[(d["MAF"] >= 0.01) & (d["INFO"] >= 0.90)].copy()

def rename_for_ldsc(d):
    m = {
        "CHROM": "CHR",
        "POS": "POS",
        "ID": "SNP",
        "A1": "A1",
        "A2": "A2",
        "BETA": "BETA",
        "SE": "SE",
        "PVAL": "P",
    }
    return d.rename(columns={k: v for k, v in m.items() if k in d.columns})

def dropna_required(d):
    req = ["SNP","A1","A2","FRQ","N","BETA","SE","P","CHR","POS","INFO"]
    return d.dropna(subset=[c for c in req if c in d.columns]).copy()

def write_ldsc(d, p):
    cols = ["SNP","A1","A2","FRQ","N","BETA","SE","P","CHR","POS","INFO"]
    p.parent.mkdir(parents=True, exist_ok=True)
    d[cols].to_csv(p, sep="\t", index=False)

def count(lab, d):
    print(f"{lab}: {len(d):,}")

def parse_args():
    ap = argparse.ArgumentParser(description="QC + LDSC format for SCZ (PGC3 wave3)")
    ap.add_argument("--in",  dest="inp",  default="Data/SCZ/PGC3_SCZ_wave3.cleaned.tsv")
    ap.add_argument("--out", dest="out", default="Data/SCZ/post-qc/PGC3_SCZ_wave3.cleaned_ldsc_ready.tsv")
    return ap.parse_args()

if __name__ == "__main__":
    a = parse_args()
    src = Path(a.inp)
    out = Path(a.out)
    d = load_tsv(src)
    count("loaded", d)
    d = to_num(d, ["CHROM","POS","BETA","SE","PVAL","IMPINFO","NCAS","NCON","FCAS","FCON"])
    d = exclude_mhc(d, "CHROM", "POS")
    count("after_exclude_MHC", d)
    d = keep_snps(d, "A1", "A2")
    count("after_remove_indels", d)
    d = drop_palindromes(d, "A1", "A2")
    count("after_remove_palindromes", d)
    d = compute_eaf(d)
    d = add_info_and_n(d)
    d = filter_maf_info(d)
    count("after_maf_info", d)
    d = rename_for_ldsc(d)
    d = dropna_required(d)
    count("after_dropna_required", d)
    write_ldsc(d, out)
    print("wrote:", out)
