#!/usr/bin/env python3
import pandas as pd
import os
import sys

cfdr_path = sys.argv[1]

ad_path = "/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv"
scz_path = "/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv"

cfdr = pd.read_csv(cfdr_path, sep="\t")
ad = pd.read_csv(ad_path, sep="\t")[["SNP", "P"]].rename(columns={"P": "P_AD"})
scz = pd.read_csv(scz_path, sep="\t")[["SNP", "P"]].rename(columns={"P": "P_SCZ"})
merged = cfdr.merge(ad, on="SNP", how="left").merge(scz, on="SNP", how="left")

base_dir = os.path.dirname(cfdr_path)
out_dir = os.path.join(base_dir, "pval_adjusted")
os.makedirs(out_dir, exist_ok=True)

out_name = os.path.basename(cfdr_path).replace(".tsv", "_pvals.tsv")
merged.to_csv(os.path.join(out_dir, out_name), sep="\t", index=False)


# python /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/scr/LD4Hits/map_original_pval.py AD_SCZ_cfdr_results_mapped.tsv
# python /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/scr/LD4Hits/map_original_pval.py AD_LONG_cfdr_results_mapped.tsv
# python /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/scr/LD4Hits/map_original_pval.py SCZ_LONG_cfdr_results_mapped.tsv

# python /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/scr/LD4Hits/map_original_pval.py AD_SCZ_shared_conjFDR0.05.tsv
# python /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/scr/LD4Hits/map_original_pval.py AD_LONG_SCZ_LONG_shared_conjFDR0.05.tsv

