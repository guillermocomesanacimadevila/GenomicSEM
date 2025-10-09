#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import sys

BASE = Path("/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data")
SZ_LDSC = BASE / "SZ" / "SZ.ldsc.sumstats"
AD_LDSC = BASE / "AD" / "AD.ldsc.sumstats"
SZ_MAP = BASE / "SZ" / "PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv"
AD_MAP = BASE / "AD" / "Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv"
OUT_SZ = BASE / "SZ" / "SZ.noMHC.ldsc.sumstats"
OUT_AD = BASE / "AD" / "AD.noMHC.ldsc.sumstats"
MHC_CHR, MHC_START, MHC_END = 6, 25_000_000, 34_000_000

def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str)

def to_num(s):
    return pd.to_numeric(s, errors="coerce")

def build_map(map_path):
    df = read_tsv(map_path)
    need = {"SNP", "CHR", "BP"}
    if not need.issubset(df.columns):
        raise ValueError(f"{map_path.name} missing columns {need - set(df.columns)}")
    m = df[["SNP", "CHR", "BP"]].dropna(subset=["SNP"]).copy()
    m["CHR"] = to_num(m["CHR"]).astype("Int64")
    m["BP"] = to_num(m["BP"]).astype("Int64")
    m = m.dropna(subset=["CHR", "BP"])
    m = m.drop_duplicates(subset=["SNP"], keep="first")
    return m

def load_ldsc_sumstats(path, require_info=False):
    df = read_tsv(path)
    need = {"SNP", "A1", "A2", "BETA", "SE", "P", "MAF"}
    if require_info:
        need = need | {"INFO"}
    if not need.issubset(df.columns):
        raise ValueError(f"{path.name} missing columns {need - set(df.columns)}")
    return df

def drop_mhc(df_with_pos):
    chrn = to_num(df_with_pos["CHR"])
    bpn = to_num(df_with_pos["BP"])
    keep = ~((chrn == MHC_CHR) & (bpn.between(MHC_START, MHC_END)))
    return df_with_pos.loc[keep]

def filter_one(label, ldsc_path, map_path, out_path, require_info=False):
    ss = load_ldsc_sumstats(ldsc_path, require_info=require_info)
    total = len(ss)
    mp = build_map(map_path)
    merged = ss.merge(mp, on="SNP", how="left", validate="m:1")
    with_pos = merged.dropna(subset=["CHR", "BP"])
    n_missing_pos = total - len(with_pos)
    after_mhc = drop_mhc(with_pos)
    n_removed_mhc = len(with_pos) - len(after_mhc)

    cols = ["SNP", "A1", "A2", "BETA", "SE", "P", "MAF"]
    if require_info and "INFO" in after_mhc.columns:
        cols.append("INFO")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    after_mhc.loc[:, cols].to_csv(out_path, sep="\t", index=False)

    print(f"[{label}] total: {total:,}")
    print(f"[{label}] without position (dropped): {n_missing_pos:,}")
    print(f"[{label}] removed in MHC (chr{MHC_CHR}:{MHC_START:,}-{MHC_END:,}): {n_removed_mhc:,}")
    print(f"[{label}] written: {len(after_mhc):,} → {out_path}")

def main():
    for p in [SZ_LDSC, AD_LDSC, SZ_MAP, AD_MAP]:
        if not p.exists():
            print(f"Missing: {p}", file=sys.stderr)
            sys.exit(1)
    filter_one("SCZ", SZ_LDSC, SZ_MAP, OUT_SZ, require_info=True)
    filter_one("AD",  AD_LDSC, AD_MAP, OUT_AD, require_info=False)

if __name__ == "__main__":
    main()
