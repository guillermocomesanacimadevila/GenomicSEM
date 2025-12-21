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

def norm_base(s):
    return s.astype(str).str.strip().str.upper()

def add_a1_a2(d):
    a1 = norm_base(d["allele1"])
    ref = norm_base(d["ref"])
    alt = norm_base(d["alt"])
    a2 = np.where(a1 == ref, alt, np.where(a1 == alt, ref, np.nan))
    d["A1"] = a1
    d["A2"] = a2
    return d

def keep_snps(d, a1c="A1", a2c="A2"):
    b = {"A", "C", "G", "T"}
    a1 = norm_base(d[a1c])
    a2 = norm_base(d[a2c])
    ok = a1.isin(b) & a2.isin(b) & (a1.str.len() == 1) & (a2.str.len() == 1)
    return d[ok].copy()

def drop_palindromes(d, a1c="A1", a2c="A2"):
    a1 = norm_base(d[a1c])
    a2 = norm_base(d[a2c])
    pal = ((a1 == "A") & (a2 == "T")) | ((a1 == "T") & (a2 == "A")) | ((a1 == "C") & (a2 == "G")) | ((a1 == "G") & (a2 == "C"))
    return d[~pal].copy()

def exclude_mhc(d, chr_col="chr", pos_col="hg38"):
    c = pd.to_numeric(d[chr_col], errors="coerce")
    p = pd.to_numeric(d[pos_col], errors="coerce")
    mhc = (c == 6) & (p >= 25_000_000) & (p <= 34_000_000)
    return d[~mhc].copy()

def format_ldsc(d, pos_build="hg38"):
    out = pd.DataFrame({
        "SNP": d["rs"],
        "A1": norm_base(d["A1"]),
        "A2": norm_base(d["A2"]),
        "FRQ": np.nan,
        "N": np.nan,
        "BETA": pd.to_numeric(d["beta"], errors="coerce"),
        "SE": pd.to_numeric(d["se"], errors="coerce"),
        "P": pd.to_numeric(d["pvalue"], errors="coerce"),
        "CHR": pd.to_numeric(d["chr"], errors="coerce"),
        "POS": pd.to_numeric(d[pos_build], errors="coerce"),
    })
    return out

def dropna_required(d):
    req = ["SNP","A1","A2","BETA","SE","P","CHR","POS"]
    return d.dropna(subset=req).copy()

def write_tsv(d, p):
    p.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(p, sep="\t", index=False)

def count(lab, d):
    print(f"{lab}: {len(d):,}")

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--pos-build", dest="pos_build", choices=["hg19","hg38"], default="hg38")
    ap.add_argument("--keep-mhc", action="store_true")
    return ap.parse_args()

if __name__ == "__main__":
    a = parse_args()
    src = Path(a.inp)
    out = Path(a.out)
    d = load_tsv(src)
    count("loaded", d)
    d = to_num(d, ["chr","hg19","hg38","beta","se","pvalue"])
    if not a.keep_mhc:
        d = exclude_mhc(d, "chr", a.pos_build)
        count("after_exclude_MHC", d)
    d = add_a1_a2(d)
    d = keep_snps(d, "A1", "A2")
    count("after_keep_snps", d)
    d = drop_palindromes(d, "A1", "A2")
    count("after_drop_palindromes", d)
    o = format_ldsc(d, a.pos_build)
    o = dropna_required(o)
    count("after_dropna_required", o)
    write_tsv(o, out)
    print("wrote:", out)
