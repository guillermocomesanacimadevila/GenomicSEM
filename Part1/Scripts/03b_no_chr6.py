#!/usr/bin/env python3
import pandas as pd
import os

VALID_ALLELES = {"A", "C", "G", "T"}
PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

def _read(path):
    return pd.read_csv(path, sep="\t", engine="python")

def _save(df, path_out):
    df.to_csv(path_out, sep="\t", index=False)
    print(f"Saved: {path_out} (n={len(df):,})")

def _rename_cols(df):
    rename_map = {
        "Effect": "A1", "Non_Effect": "A2",
        "EA": "A1", "NEA": "A2",
        "ALT": "A1", "REF": "A2",
        "Effect_allele_freq": "MAF", "Freq1": "MAF", "FREQ": "MAF",
        "ImpInfo": "INFO", "IMPINFO": "INFO",
        "Beta": "BETA", "beta": "BETA", "Effect_size": "BETA", "OR": "BETA",
        "StdErr": "SE", "SE": "SE",
        "PVAL": "P", "Pvalue": "P", "p": "P",
        "ID": "SNP"
    }
    return df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

def _basic_qc(df, info_cut=0.90, maf_cut=0.01):
    if "A1" not in df.columns or "A2" not in df.columns:
        raise ValueError("Missing allele columns (A1/A2) after renaming.")
    df["A1"] = df["A1"].astype(str).str.upper()
    df["A2"] = df["A2"].astype(str).str.upper()
    df = df[df["A1"].isin(VALID_ALLELES) & df["A2"].isin(VALID_ALLELES)].copy()
    pal = df.apply(lambda r: (r["A1"], r["A2"]) in PALINDROMIC, axis=1)
    if "MAF" in df.columns:
        df["MAF"] = pd.to_numeric(df["MAF"], errors="coerce")
        ambiguous = pal & df["MAF"].between(0.4, 0.6, inclusive="both")
        df = df[~ambiguous].copy()
    else:
        df = df[~pal].copy()
    if "INFO" in df.columns:
        df["INFO"] = pd.to_numeric(df["INFO"], errors="coerce")
        df = df[df["INFO"] >= info_cut]
    if "MAF" in df.columns:
        df = df[df["MAF"] >= maf_cut]
    required_any = [c for c in ["SNP", "BETA", "SE", "P"] if c in df.columns]
    if not required_any:
        raise ValueError("Missing core GWAS columns (SNP/BETA/SE/P) after renaming.")
    df = df.dropna(subset=list(set(["SNP", "A1", "A2"] + required_any)))
    return df

def drop_chr6_from_prepared(prepared_path, out_path=None, chr_col="CHR"):
    df_raw = _read(prepared_path)
    df_raw = _rename_cols(df_raw)
    if chr_col not in df_raw.columns:
        raise ValueError(f"{prepared_path}: missing '{chr_col}' column.")

    before_total = len(df_raw)
    chr_raw = pd.to_numeric(df_raw[chr_col], errors="coerce")
    before_chr6 = int((chr_raw == 6).sum())

    df_qc = _basic_qc(df_raw.copy(), info_cut=0.90, maf_cut=0.01)
    after_qc = len(df_qc)

    df_qc[chr_col] = pd.to_numeric(df_qc[chr_col], errors="coerce")
    df_out = df_qc[df_qc[chr_col] != 6].copy()
    after_chr6 = len(df_out)

    if out_path is None:
        base, ext = os.path.splitext(prepared_path)
        ext = ext if ext else ".tsv"
        out_path = f"{base}.clean.noChr6{ext}"

    _save(df_out, out_path)

    print(f"RAW total SNPs: {before_total:,}")
    print(f"RAW chr6 SNPs (before QC): {before_chr6:,}")
    print(f"SNPs after QC: {after_qc:,}")
    print(f"Dropped on chr6 (post-QC): {after_qc - after_chr6:,}")
    print(f"FINAL SNP COUNT (QC + no chr6): {after_chr6:,}\n")
    return out_path

def drop_chr6_from_ldsc(ldsc_path, prepared_path, out_path=None, chr_col="CHR", snp_col="SNP"):
    ldsc_raw = _read(ldsc_path)
    prep_raw = _read(prepared_path)
    ldsc_raw = _rename_cols(ldsc_raw)
    prep_raw = _rename_cols(prep_raw)
    if snp_col not in ldsc_raw.columns:
        raise ValueError(f"{ldsc_path}: missing '{snp_col}'")
    if chr_col not in prep_raw.columns:
        raise ValueError(f"{prepared_path}: missing '{chr_col}'")

    before_total = len(ldsc_raw)
    prep_raw[chr_col] = pd.to_numeric(prep_raw[chr_col], errors="coerce")
    chr6_rsids_raw = set(prep_raw.loc[prep_raw[chr_col] == 6, snp_col].astype(str))
    ldsc_raw[snp_col] = ldsc_raw[snp_col].astype(str)
    before_chr6 = int(ldsc_raw[snp_col].isin(chr6_rsids_raw).sum())

    ldsc_qc = _basic_qc(ldsc_raw.copy(), info_cut=0.90, maf_cut=0.01)
    after_qc = len(ldsc_qc)

    # recompute chr6 rsids from the same prepared file (position source of truth)
    chr6_rsids = chr6_rsids_raw
    ldsc_out = ldsc_qc[~ldsc_qc[snp_col].isin(chr6_rsids)].copy()
    after_chr6 = len(ldsc_out)

    if out_path is None:
        base, ext = os.path.splitext(ldsc_path)
        ext = ext if ext else ".sumstats"
        out_path = f"{base}.clean.noChr6{ext}"

    _save(ldsc_out, out_path)

    print(f"RAW total SNPs: {before_total:,}")
    print(f"RAW chr6 SNPs in LDSC file (via rsID ∩ chr6 from prepared): {before_chr6:,}")
    print(f"SNPs after QC: {after_qc:,}")
    print(f"Dropped on chr6 (post-QC): {after_qc - after_chr6:,}")
    print(f"FINAL SNP COUNT (QC + no chr6): {after_chr6:,}\n")
    return out_path

if __name__ == "__main__":
    sz_ldsc = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/SZ/SZ.ldsc.sumstats"
    ad_ldsc = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/AD/AD.ldsc.sumstats"
    sz_prepared = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv"
    ad_prepared = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv"

    print("=== Prepared-for-plots → drop chr6 with full QC ===")
    sz_p_nochr6 = drop_chr6_from_prepared(sz_prepared)
    ad_p_nochr6 = drop_chr6_from_prepared(ad_prepared)

    print("=== LDSC inputs → drop chr6 with full QC ===")
    sz_ldsc_nochr6 = drop_chr6_from_ldsc(sz_ldsc, sz_prepared)
    ad_ldsc_nochr6 = drop_chr6_from_ldsc(ad_ldsc, ad_prepared)

    print("Use these in R:")
    print("  SZ noChr6 (LDSC):", sz_ldsc_nochr6)
    print("  AD noChr6 (LDSC):", ad_ldsc_nochr6)
