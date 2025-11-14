#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import os
import sys

PATH = Path("/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs")
REF_DIR = Path(os.environ.get("REF_DIR", "/Users/c24102394/ref/ldsc/1000G_EUR_Phase3_plink"))
REF_PATTERN = os.environ.get("REF_PATTERN", "1000G.EUR.QC.*.bim")

SHOW_LDSC = True
SHOW_LAVA = True
SHOW_CONJFDR = True

P_ASSET_THRESHOLD = 0.05

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 160)
pd.set_option("display.colheader_justify", "center")

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR / "conjFDR"))
from map_hits import read_ref, map_one_hits

def print_section(title):
    print("\n" + "=" * 70)
    print(title.upper())
    print("=" * 70)

def clean_ldsc_table(path):
    if not path.exists():
        print(f"Missing file: {path}")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", comment="#")
    if df.shape[1] == 1:
        df = pd.read_csv(path, sep=",", comment="#")
    df.columns = df.columns.str.replace('"', '').str.strip()
    df = df.rename(columns={c: c.strip() for c in df.columns})
    if "Unnamed: 0" in df.columns:
        df = df.drop(columns=["Unnamed: 0"])
    return df

def load_ldsc_results():
    p = PATH / "ldsc"
    files = {
        "AD-SCZ": p / "ad_scz/ad_scz_ldsc_results.tsv",
        "AD-LONG": p / "ad_age/ad_age_ldsc_results.tsv",
        "SCZ-LONG": p / "scz_age/scz_age_ldsc_results.tsv",
    }
    out = {}
    for name, f in files.items():
        out[name] = clean_ldsc_table(f)
    return out

def show_ldsc_outputs():
    res = load_ldsc_results()
    for label, df in res.items():
        print_section(f"LDSC genome-wide results for {label}")
        if df.empty:
            print("No data found or unreadable format.")
            continue
        print("\nResults:")
        print(df.to_string(index=False))
        print("\n")

def show_lava_hits():
    p = PATH / "lava" / "ad_scz_age"
    lava_file = p / "LAVA_local_rg_bivariate.tsv"
    if not lava_file.exists():
        print_section("LAVA local genetic correlation hits (bivariate)")
        print("Missing LAVA file:", lava_file)
        return
    df = pd.read_csv(lava_file, sep="\t")
    print_section("LAVA local genetic correlation hits (bivariate)")
    def one(pair_code, label):
        sub = df[df["pair"] == pair_code]
        if sub.empty:
            print(f"{label}: no regions")
            return
        n_raw = (sub["p"] < 0.05).sum() if "p" in sub.columns else 0
        n_fdr = (sub["q_fdr"] < 0.05).sum() if "q_fdr" in sub.columns else 0
        n_bonf = (sub["p_bonf_paper"] < 0.05).sum() if "p_bonf_paper" in sub.columns else 0
        print(f"\n{label}:")
        print(f"  raw p<0.05 = {n_raw}")
        print(f"  FDR<0.05  = {n_fdr}")
        print(f"  Bonf<0.05 = {n_bonf}")
        bonf_hits = sub[sub.get("p_bonf_paper", 1) < 0.05] if "p_bonf_paper" in sub.columns else pd.DataFrame()
        if not bonf_hits.empty:
            print("\n  Bonferroni-significant loci:")
            print(bonf_hits.to_string(index=False))
            print("\n")
    one("AD_SCZ", "AD-SCZ")
    one("AD_AGE", "AD-LONG")
    one("SCZ_AGE", "SCZ-LONG")

def safe_read_ref_recursive(ref_dir, pattern):
    ref_dir = Path(ref_dir)
    bim_files = sorted(ref_dir.glob(pattern))
    if len(bim_files) == 0:
        bim_files = sorted(ref_dir.rglob("*.bim"))
    if len(bim_files) == 0:
        print_section("Reference mapping")
        print(f"No .bim files found in {ref_dir} (pattern {pattern} or recursive *.bim)")
        return pd.DataFrame(columns=["SNP", "CHR", "BP"])
    return read_ref(ref_dir)

def map_and_print(df, ref, mapped_dir, label):
    if df is None or df.empty:
        print("None")
        return
    mapped_dir.mkdir(parents=True, exist_ok=True)
    tmp = mapped_dir / f"{label}.tsv"
    df.to_csv(tmp, sep="\t", index=False)
    if ref.empty:
        print(df.to_string(index=False))
        return
    m = map_one_hits(tmp, ref)
    cols = [c for c in ["SNP","CHR","BP","START","END"] if c in m.columns] + [c for c in m.columns if c not in ["SNP","CHR","BP","START","END"]]
    print(m[cols].to_string(index=False))

def show_conjfdr_results():
    p = PATH / "conjFDR"
    mapped_dir = p / "mapped"
    print_section("conjFDR results per trait pair")
    pairs = {
        "AD-SCZ": "AD_SCZ",
        "AD-LONG": "AD_LONG",
        "SCZ-LONG": "SCZ_LONG",
    }
    ad_scz_cf = None
    ad_scz_hits = None
    for label, prefix in pairs.items():
        cf_path = p / f"{prefix}_cfdr_results.tsv"
        hits_path = p / f"{prefix}_shared_hits.tsv"
        if not cf_path.exists() or not hits_path.exists():
            continue
        cf = pd.read_csv(cf_path, sep="\t")
        hits = pd.read_csv(hits_path, sep="\t")
        n_total = cf.shape[0]
        n_conj05 = (cf["conj_fdr"] < 0.05).sum() if "conj_fdr" in cf.columns else 0
        n_hits = hits.shape[0]
        print(f"\n{label}:")
        print(f"  total SNPs       = {n_total}")
        print(f"  conj_fdr < 0.05  = {n_conj05}")
        print(f"  shared_hits file = {n_hits}")
        if label == "AD-SCZ":
            ad_scz_cf = cf
            ad_scz_hits = hits
    ref = safe_read_ref_recursive(REF_DIR, REF_PATTERN)
    if ad_scz_cf is not None:
        print_section("AD-SCZ SNPs with conj_fdr < 0.05 (mapped to reference)")
        cf_top = ad_scz_cf[ad_scz_cf["conj_fdr"] < 0.05].copy()
        if not cf_top.empty:
            cf_top["SNP"] = cf_top["SNP"].astype(str)
            map_and_print(cf_top, ref, mapped_dir, "AD_SCZ_cfdr_lt_0.05")
        else:
            print("None")
    if ad_scz_hits is not None:
        print_section("AD-SCZ shared_hits table mapped to reference")
        if not ad_scz_hits.empty:
            ad_scz_hits["SNP"] = ad_scz_hits["SNP"].astype(str)
            map_and_print(ad_scz_hits, ref, mapped_dir, "AD_SCZ_shared_hits")
        else:
            print("None")
    print_section("Shared SNPs across conjFDR trait pairs")
    shared_path = mapped_dir / "shared_snps_across_trait_pairs_full.csv"
    if shared_path.exists():
        shared = pd.read_csv(shared_path)
        n_shared = shared.shape[0]
        print(f"SNPs present in ≥2 trait pairs: {n_shared}\n")
        if not shared.empty:
            print("\nShared SNP table:\n")
            print(shared.to_string(index=False))
            print("\n")
    else:
        print("No shared_snps_across_trait_pairs_full.csv found")

def show_asset():
    p = PATH / "ASSET"
    ad_scz_age = p / "ad_scz_age"
    print_section("ASSET")
    if not ad_scz_age.exists():
        print("ASSET directory not found")
        return
    for f in os.listdir(ad_scz_age):
        if f.endswith(".tsv"):
            print(f"\n{f}:")
    pleio_path = ad_scz_age / "AD_SCZ_AGE_ASSET_pleiotropic_hits.tsv"
    full_path = ad_scz_age / "AD_SCZ_AGE_ASSET_results.tsv"
    pleio = pd.read_csv(pleio_path, sep="\t") if pleio_path.exists() else pd.DataFrame()
    full = pd.read_csv(full_path, sep="\t") if full_path.exists() else pd.DataFrame()
    n_total = full.shape[0]
    n_asset_pleio_std = (pleio["P_ASSET"] < 0.05).sum() if "P_ASSET" in pleio.columns else 0
    n_asset_full_std = (full["P_ASSET"] < 0.05).sum() if "P_ASSET" in full.columns else 0
    n_asset_full_relaxed = (full["P_ASSET"] < P_ASSET_THRESHOLD).sum() if "P_ASSET" in full.columns else 0
    print(f"\nAsset sumary:")
    print(f"Total SNPs = {n_total}")
    print(f"Number of Pleiotropic SNPs (all traits, P_ASSET<0.05) = {n_asset_pleio_std}")
    print(f"Total SNPs below P_ASSET<0.05 = {n_asset_full_std}")
    print(f"Total SNPs below P_ASSET<{P_ASSET_THRESHOLD} = {n_asset_full_relaxed}")
    if "P_ASSET" in full.columns and n_asset_full_relaxed > 0:
        sig = full[full["P_ASSET"] < P_ASSET_THRESHOLD]
        print("\nSNPs with P_ASSET below relaxed threshold:\n")
        print(sig.to_string(index=False))

def main():
    if SHOW_LDSC:
        show_ldsc_outputs()
    if SHOW_LAVA:
        show_lava_hits()
    if SHOW_CONJFDR:
        show_conjfdr_results()

if __name__ == "__main__":
    main()
