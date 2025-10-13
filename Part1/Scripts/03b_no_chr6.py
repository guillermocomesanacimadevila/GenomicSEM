import pandas as pd
import os

def prepared_to_ldsc_nochr6(path_in, out_path=None, drop_chr=6):
    df = pd.read_csv(path_in, sep="\t", engine="python")
    colmap = {c: c.strip() for c in df.columns}
    df.rename(columns=colmap, inplace=True)
    required_in = ["CHR", "SNP", "Effect", "Non_Effect", "Beta", "SE", "P"]
    missing = [c for c in required_in if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path_in}: {missing}")

    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df = df[df["CHR"] != drop_chr].copy()
    out = pd.DataFrame({
        "SNP":  df["SNP"],
        "A1":   df["Effect"],
        "A2":   df["Non_Effect"],
        "BETA": pd.to_numeric(df["Beta"], errors="coerce"),
        "SE":   pd.to_numeric(df["SE"], errors="coerce"),
        "P":    pd.to_numeric(df["P"], errors="coerce"),
    })

    if "IMPINFO" in df.columns:
        out["INFO"] = pd.to_numeric(df["IMPINFO"], errors="coerce")
    else:
        out["INFO"] = pd.NA

    if "maf" in df.columns:
        out["MAF"] = pd.to_numeric(df["maf"], errors="coerce")
    else:
        eaf_col = None
        for cand in ("eaf", "Effect_allele_freq"):
            if cand in df.columns:
                eaf_col = cand
                break
        if eaf_col is not None:
            eaf = pd.to_numeric(df[eaf_col], errors="coerce")
            out["MAF"] = eaf.where(eaf <= 0.5, 1 - eaf)
        else:
            out["MAF"] = pd.NA

    out = out[["SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "INFO"]]
    out = out.dropna(subset=["SNP", "A1", "A2", "BETA", "SE", "P"])
    if out_path is None:
        base, _ext = os.path.splitext(path_in)
        out_path = f"{base}_no_chr6.ldsc.sumstats"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Saved LDSC-style sumstats (CHR≠6) to: {out_path}")
    return out_path

if __name__ == "__main__":
    # path = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv"
    path = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv"
    prepared_to_ldsc_nochr6(path)
