import pandas as pd
import os

def drop_chr6_from_ldsc(ldsc_path, prepared_path, out_path=None):
    ldsc = pd.read_csv(ldsc_path, sep="\t", engine="python")
    prep = pd.read_csv(prepared_path, sep="\t", engine="python")

    for need in [["SNP","MAF"], ["SNP","CHR"]]:
        src = ldsc if "MAF" in need else prep
        missing = [c for c in need if c not in src.columns]
        if missing:
            raise ValueError(f"Missing columns {missing} in {'LDSC' if src is ldsc else 'prepared'} file")

    prep["CHR"] = pd.to_numeric(prep["CHR"], errors="coerce")
    chr6 = set(prep.loc[prep["CHR"] == 6, "SNP"].astype(str))
    ldsc["SNP"] = ldsc["SNP"].astype(str)
    out = ldsc[~ldsc["SNP"].isin(chr6)].copy()
    if out_path is None:
        base, ext = os.path.splitext(ldsc_path)
        out_path = f"{base}.noChr6{ext}"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Saved: {out_path}  (kept {len(out):,} SNPs)")
    return out_path

if __name__ == "__main__":
    sz_ldsc = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/SZ.ldsc.sumstats"
    ad_ldsc = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/AD.ldsc.sumstats"
    sz_prepared = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv"
    ad_prepared = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv"

    drop_chr6_from_ldsc(sz_ldsc, sz_prepared)  # -> SZ/SZ.ldsc.noChr6.sumstats
    drop_chr6_from_ldsc(ad_ldsc, ad_prepared)  # -> AD/AD.ldsc.noChr6.sumstats
