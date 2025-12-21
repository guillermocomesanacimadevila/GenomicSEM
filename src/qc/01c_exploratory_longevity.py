#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

def load_tsv(p):
    return pd.read_csv(p, sep="\t", dtype=str, low_memory=False)

def rename_cols(d):
    d = d.drop(
        columns=[
            "snpid",
            "beta1_healthspan","se_healthspan",
            "beta1_lifespan","se_lifespan",
            "beta1_longevity","se_longevity",
        ],
        errors="ignore"
    )
    return d.reset_index(drop=True).rename(
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

def to_num(d):
    for c in ["CHR","POS","N","BETA","SE","P","EAF","INFO"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d

def exclude_mhc(d):
    mhc = (d["CHR"] == 6) & d["POS"].between(25_000_000, 34_000_000)
    return d[~mhc].copy()

def keep_snps(d):
    b = {"A","C","G","T"}
    a1 = d["A1"].astype(str).str.upper()
    a2 = d["A2"].astype(str).str.upper()
    snp = a1.isin(b) & a2.isin(b) & (a1.str.len() == 1) & (a2.str.len() == 1)
    gap = a1.str.contains("-", na=False) | a2.str.contains("-", na=False)
    return d[snp & ~gap].copy()

def drop_palindromes(d):
    a1 = d["A1"].astype(str).str.upper()
    a2 = d["A2"].astype(str).str.upper()
    pal = ((a1 == "A") & (a2 == "T")) | ((a1 == "T") & (a2 == "A")) | ((a1 == "C") & (a2 == "G")) | ((a1 == "G") & (a2 == "C"))
    return d[~pal].copy()

def maf_info_filter(d):
    d["FRQ"] = d["EAF"]
    d["MAF"] = np.minimum(d["FRQ"], 1 - d["FRQ"])
    return d[(d["MAF"] >= 0.01) & (d["INFO"] >= 0.90)].copy()

def dropna_req(d):
    req = ["SNP","A1","A2","FRQ","N","BETA","SE","P","CHR","POS","INFO"]
    return d.dropna(subset=[c for c in req if c in d.columns]).copy()

def write_ldsc(d, p):
    cols = ["SNP","A1","A2","FRQ","N","BETA","SE","P","CHR","POS","INFO"]
    p.parent.mkdir(parents=True, exist_ok=True)
    d[cols].to_csv(p, sep="\t", index=False)

def count(lab, d):
    print(f"{lab}: {len(d):,}")

def parse_args():
    ap = argparse.ArgumentParser(description="QC + LDSC format for Timmers 2020 healthspan/lifespan/longevity")
    ap.add_argument("--in", dest="inp", default="Data/AGE/timmers2020_healthspan_lifespan_longevity.tsv")
    ap.add_argument("--out", dest="out", default="Data/AGE/post-qc/AGE_ldsc_ready.tsv")
    return ap.parse_args()

if __name__ == "__main__":
    a = parse_args()
    src = Path(a.inp)
    out = Path(a.out)
    d = load_tsv(src)
    count("loaded", d)
    d = rename_cols(d)
    d = to_num(d)
    count("after_rename_numeric", d)
    d = exclude_mhc(d)
    count("after_exclude_MHC", d)
    d = keep_snps(d)
    count("after_remove_indels", d)
    d = drop_palindromes(d)
    count("after_remove_palindromes", d)
    d = maf_info_filter(d)
    count("after_maf_info", d)
    d = dropna_req(d)
    count("after_dropna_required", d)
    write_ldsc(d, out)
    print("wrote:", out)