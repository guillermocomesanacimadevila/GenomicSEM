#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd

def convert(path: str) -> str:
    df = pd.read_csv(path, sep=r"\s+")
    tsv_path = os.path.splitext(path)[0] + ".tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path

def qc_checks(df: pd.DataFrame, drop_indels: bool = False) -> pd.DataFrame:
    print("===== QC Summary =====")
    print(f"Number of SNPs: {df.shape[0]:,}")
    print(f"Number of columns: {df.shape[1]}")
    print("\nColumn names:", list(df.columns))
    print("\nMissing values per column:")
    print(df.isna().sum())

    cols_lower = [c.lower() for c in df.columns]
    if "beta" in cols_lower:
        print("\nEffect size column: BETA")
    elif "or" in cols_lower:
        print("\nEffect size column: OR")
    else:
        print("\nNo BETA/OR column found")

    pcol = "P" if "P" in df.columns else None
    if pcol:
        vmax = pd.to_numeric(df[pcol], errors="coerce").max()
        print(f"{pcol} column looks like {'-log10(P)' if (pd.notna(vmax) and vmax > 1) else 'raw P-values'}")
    else:
        print("No P/PVAL column found")

    allele_pairs = [
        ("Effect", "Non_Effect"),
        ("A1", "A2"),
        ("EA", "NEA"),
        ("ALT", "REF"),
        ("ALLELE1", "ALLELE2"),
    ]
    a1, a2 = None, None
    for c1, c2 in allele_pairs:
        if c1 in df.columns and c2 in df.columns:
            a1, a2 = c1, c2
            break

    if a1 and a2:
        bases = {"A", "C", "G", "T"}
        a1s = df[a1].astype(str).str.upper()
        a2s = df[a2].astype(str).str.upper()
        is_single = a1s.isin(bases) & a2s.isin(bases) & (a1s.str.len() == 1) & (a2s.str.len() == 1)
        has_gap = a1s.str.contains("-", na=False) | a2s.str.contains("-", na=False)
        is_snp = is_single & (~has_gap)
        is_indel = ~is_snp
        n_indels = int(is_indel.sum())
        n_snps_before = int(is_snp.sum())
        print(f"\nINDELs detected: {n_indels:,}")
        print(f"SNPs (before drop): {n_snps_before:,}")
        if drop_indels:
            df = df[is_snp].copy()
            print(f"SNPs (after drop): {df.shape[0]:,}")
    else:
        print("\nINDEL check skipped: no allele columns found")

    print("=======================")
    return df

def prepare_for_plots(df: pd.DataFrame):
    df_out = df.copy()

    def _norm_chr(x):
        if pd.isna(x): return np.nan
        s = str(x).strip()
        s = re.sub(r"^chr", "", s, flags=re.IGNORECASE).upper()
        if s == "X": return 23
        if s == "Y": return 24
        if s in {"MT", "M", "MITO"}: return 25
        try:
            return int(float(s))
        except Exception:
            return np.nan

    df_out["CHR_int"] = df_out["CHR"].map(_norm_chr).astype("Int64")
    df_out["BP_int"] = pd.to_numeric(df_out["BP"], errors="coerce").astype("Int64")

    pvals = pd.to_numeric(df_out["P"], errors="coerce")
    tiny = 1e-300
    pvals = pvals.where((pvals > 0) & (pvals <= 1), np.nan).clip(lower=tiny)
    df_out["minus_log10p"] = -np.log10(pvals)

    eff = df_out["Effect"].astype(str)
    non = df_out["Non_Effect"].astype(str)
    bases = {"A", "C", "G", "T"}

    is_snp = eff.isin(bases) & non.isin(bases) & (eff.str.len() == 1) & (non.str.len() == 1)
    is_indel = ~is_snp | eff.str.contains("-", na=False) | non.str.contains("-", na=False)

    ambiguous = (is_snp & (
        ((eff == "A") & (non == "T")) |
        ((eff == "T") & (non == "A")) |
        ((eff == "C") & (non == "G")) |
        ((eff == "G") & (non == "C"))
    ))

    df_out["is_snp"] = is_snp.fillna(False)
    df_out["is_indel"] = is_indel.fillna(False)
    df_out["is_ambiguous"] = ambiguous.fillna(False)

    snp_id = df_out["SNP"].astype(str)
    rs_mask = snp_id.str.startswith("rs", na=False)
    fallback = (
        df_out["CHR_int"].astype(str)
        + ":" + df_out["BP_int"].astype(str)
        + ":" + non + "_" + eff
    )
    df_out["variant_id"] = np.where(
        rs_mask, snp_id,
        np.where(df_out["Test_MarkerName"].notna(),
                 df_out["Test_MarkerName"].astype(str),
                 fallback)
    )

    if "Effect_allele_freq" in df_out.columns:
        eaf = pd.to_numeric(df_out["Effect_allele_freq"], errors="coerce")
        df_out["eaf"] = eaf
        df_out["maf"] = np.minimum(eaf, 1 - eaf)
    else:
        df_out["eaf"] = np.nan
        df_out["maf"] = np.nan

    before = len(df_out)
    df_out = df_out.dropna(subset=["CHR_int", "BP_int", "minus_log10p"]).copy()
    n_dropped = before - len(df_out)

    autos = sorted([c for c in df_out["CHR_int"].unique() if pd.notna(c) and c <= 22])
    others = [c for c in [23, 24, 25] if c in df_out["CHR_int"].unique()]
    chrom_order = autos + others

    df_out.sort_values(["CHR_int", "BP_int"], inplace=True)
    xticks, xlabels = [], []
    running = 0
    for c in chrom_order:
        mask = df_out["CHR_int"] == c
        if not mask.any(): continue
        bp_min, bp_max = df_out.loc[mask, "BP_int"].agg(["min", "max"])
        df_out.loc[mask, "pos_cum"] = df_out.loc[mask, "BP_int"] + running
        mid = running + (bp_max - bp_min) / 2
        xticks.append(mid)
        xlabels.append(str(c if c <= 22 else {23: "X", 24: "Y", 25: "MT"}[c]))
        running += (bp_max - bp_min + 1)

    df_out["pos_cum"] = df_out["pos_cum"].astype(np.int64)

    meta = {"n_dropped": n_dropped, "chrom_order": chrom_order}
    return df_out, xticks, xlabels, meta

if __name__ == "__main__":
    path_txt = "/Users/guillermocomesanacimadevila/Desktop/PhD/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.txt"
    df = pd.read_csv(path_txt, sep=r"\s+")
    df = qc_checks(df, drop_indels=True)

    print("\n[1/2] Writing RAW TSV (no -log10 transform)...")
    raw_tsv_path = convert(path_txt)
    print(f"RAW TSV: {raw_tsv_path}")

    print("\n[2/2] Preparing DataFrame for plots...")
    df_prep, xticks, xlabels, meta = prepare_for_plots(df)
    print("Prepared rows:", df_prep.shape[0], "| Dropped:", meta["n_dropped"])
    print("Chromosomes detected:", meta["chrom_order"])

    base_no_ext, _ = os.path.splitext(path_txt)
    prepared_tsv_path = base_no_ext + ".prepared_for_plots.tsv"
    df_prep.to_csv(prepared_tsv_path, sep="\t", index=False)
    print(f"PREPARED TSV: {prepared_tsv_path}")
    print("Done.")



# TO DO´s
# Do all of this but for the SZ GWAS
# Use Genomic SEM package for LDSC (R)
# Get genetic correlation between both GWAS - get a small table + heatmap 22
# r(theta) should be 0 or near 0 
# Also make sure Beta is not -log10(p) when slapping data into GenomicSEM LDSC