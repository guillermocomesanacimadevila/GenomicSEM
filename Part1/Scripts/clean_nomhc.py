#!/usr/bin/env python3
import pandas as pd
import os

# hg19 default
MHC_CHR = 6
MHC_START = 25_000_000
MHC_END   = 34_000_000

def _ensure_tabular_read(path):
    return pd.read_csv(path, sep="\t", engine="python")

def _save_tsv(df, path_out):
    df.to_csv(path_out, sep="\t", index=False)
    print(f"Saved: {path_out} (n={len(df):,})")


def strip_mhc_from_prepared(path_in, path_out=None,
                            chr_col="CHR", bp_col="BP",
                            chrom=MHC_CHR, start=MHC_START, end=MHC_END):
    df = _ensure_tabular_read(path_in)
    cols = list(df.columns)
    if chr_col not in df.columns or bp_col not in df.columns:
        raise ValueError(f"{path_in}: missing '{chr_col}' and/or '{bp_col}' columns.")

    df[chr_col] = pd.to_numeric(df[chr_col], errors="coerce")
    df[bp_col]  = pd.to_numeric(df[bp_col], errors="coerce")
    before = len(df)
    mask_mhc = (df[chr_col] == chrom) & (df[bp_col] >= start) & (df[bp_col] <= end)
    df_out = df.loc[~mask_mhc, :].copy()
    after = len(df_out)

    if path_out is None:
        base, ext = os.path.splitext(path_in)
        ext = ext if ext else ".tsv"
        path_out = f"{base}.noMHC{ext}"

    df_out = df_out[cols]
    _save_tsv(df_out, path_out)
    print(f"Removed {before - after:,} SNP rows in chr{chrom}:{start:,}-{end:,} (MHC).")
    return path_out

def strip_mhc_from_ldsc(ldsc_path, prepared_path, path_out=None,
                        chr_col="CHR", bp_col="BP", snp_col="SNP",
                        chrom=MHC_CHR, start=MHC_START, end=MHC_END):
    ldsc = _ensure_tabular_read(ldsc_path)
    prep = _ensure_tabular_read(prepared_path)
    if snp_col not in ldsc.columns:
        raise ValueError(f"{ldsc_path}: missing '{snp_col}' column.")
    for c in (snp_col, chr_col, bp_col):
        if c not in prep.columns:
            raise ValueError(f"{prepared_path}: missing column '{c}' required to locate MHC.")

    prep[chr_col] = pd.to_numeric(prep[chr_col], errors="coerce")
    prep[bp_col]  = pd.to_numeric(prep[bp_col], errors="coerce")
    mhc_rsids = set(prep.loc[
        (prep[chr_col] == chrom) & (prep[bp_col] >= start) & (prep[bp_col] <= end),
        snp_col
    ].astype(str))

    ldsc[snp_col] = ldsc[snp_col].astype(str)
    before = len(ldsc)
    ldsc_out = ldsc.loc[~ldsc[snp_col].isin(mhc_rsids), :].copy()
    after = len(ldsc_out)
    if path_out is None:
        base, ext = os.path.splitext(ldsc_path)
        ext = ext if ext else ".sumstats"
        path_out = f"{base}.noMHC{ext}"

    _save_tsv(ldsc_out, path_out)
    print(f"Removed {before - after:,} SNP rows by rsID match to MHC (chr{chrom}:{start:,}-{end:,}).")
    return path_out

if __name__ == "__main__":
    sz_prepared = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv"
    ad_prepared = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv"
    sz_ldsc = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/SZ.ldsc.sumstats"
    ad_ldsc = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/AD.ldsc.sumstats"

    print("=== Removing MHC from prepared_for_plots files ===")
    sz_prepared_noMHC = strip_mhc_from_prepared(sz_prepared)  # -> *.noMHC.tsv
    ad_prepared_noMHC = strip_mhc_from_prepared(ad_prepared)  # -> *.noMHC.tsv

    print("\n=== Removing MHC from LDSC files (preserving MAF/INFO) ===")
    sz_ldsc_noMHC = strip_mhc_from_ldsc(sz_ldsc, sz_prepared)  # -> *.noMHC.ldsc.sumstats
    ad_ldsc_noMHC = strip_mhc_from_ldsc(ad_ldsc, ad_prepared)  # -> *.noMHC.ldsc.sumstats

    print("\nDone. Use these in R:")
    print("  SZ noMHC (LDSC):", sz_ldsc_noMHC)
    print("  AD noMHC (LDSC):", ad_ldsc_noMHC)


# NEED to use the cleaned file - THIS ONE no-MHC
# NEED to clean WITH MHC AND ALSO NO CHR6
