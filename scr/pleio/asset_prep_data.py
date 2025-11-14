#!/usr/bin/env python3
import argparse
import pandas as pd

def read_sumstats(path, prefix):
    cols = ["SNP", "A1", "A2", "BETA", "SE", "P"]
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in cols or c.upper() in cols)
    df = df.rename(columns={c: c.upper() for c in df.columns})
    df = df[cols]
    df = df.dropna()
    df = df.rename(columns={
        "A1": f"A1_{prefix}",
        "A2": f"A2_{prefix}",
        "BETA": f"BETA_{prefix}",
        "SE": f"SE_{prefix}",
        "P": f"P_{prefix}",
    })
    return df

def harmonise_three(d1, p1, d2, p2, d3, p3):
    m = d1.merge(d2, on="SNP", how="inner").merge(d3, on="SNP", how="inner")
    a1_ref = m[f"A1_{p1}"]
    a2_ref = m[f"A2_{p1}"]

    a1_2 = m[f"A1_{p2}"]
    a2_2 = m[f"A2_{p2}"]
    same2 = (a1_2 == a1_ref) & (a2_2 == a2_ref)
    flip2 = (a1_2 == a2_ref) & (a2_2 == a1_ref)

    a1_3 = m[f"A1_{p3}"]
    a2_3 = m[f"A2_{p3}"]
    same3 = (a1_3 == a1_ref) & (a2_3 == a2_ref)
    flip3 = (a1_3 == a2_ref) & (a2_3 == a1_ref)

    keep = (same2 | flip2) & (same3 | flip3)
    m = m[keep].copy()

    flip2_keep = flip2[keep].to_numpy()
    flip3_keep = flip3[keep].to_numpy()

    m.loc[flip2_keep, f"BETA_{p2}"] = -m.loc[flip2_keep, f"BETA_{p2}"]
    m.loc[flip3_keep, f"BETA_{p3}"] = -m.loc[flip3_keep, f"BETA_{p3}"]

    m["A1"] = m[f"A1_{p1}"]
    m["A2"] = m[f"A2_{p1}"]

    out = m[[
        "SNP",
        "A1",
        "A2",
        f"BETA_{p1}",
        f"SE_{p1}",
        f"P_{p1}",
        f"BETA_{p2}",
        f"SE_{p2}",
        f"P_{p2}",
        f"BETA_{p3}",
        f"SE_{p3}",
        f"P_{p3}",
    ]].dropna()

    out = out[out["A1"] != out["A2"]]
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--t1", required=True)
    p.add_argument("--t1-prefix", required=True)
    p.add_argument("--t2", required=True)
    p.add_argument("--t2-prefix", required=True)
    p.add_argument("--t3", required=True)
    p.add_argument("--t3-prefix", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    d1 = read_sumstats(args.t1, args.t1_prefix)
    d2 = read_sumstats(args.t2, args.t2_prefix)
    d3 = read_sumstats(args.t3, args.t3_prefix)

    out = harmonise_three(d1, args.t1_prefix, d2, args.t2_prefix, d3, args.t3_prefix)
    out.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
