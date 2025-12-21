#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def calc_neff(df_path, cases, controls):
    df = pd.read_csv(df_path, sep="\t", dtype=str, low_memory=False)
    n_eff = 4 / ((1 / float(cases)) + (1 / float(controls)))
    df["N"] = round(n_eff, 0)
    pd.to_numeric(df["N"], errors="coerce")
    return df

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in",  dest="inp", required=True)
    p.add_argument("--out", dest="out", required=True)
    p.add_argument("--cases", required=True, type=float)
    p.add_argument("--controls", required=True, type=float)
    a = p.parse_args()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = calc_neff(a.inp, a.cases, a.controls)
    df.to_csv(out, sep="\t", index=False)
    print(f"written: {out}")

if __name__ == "__main__":
    main()
