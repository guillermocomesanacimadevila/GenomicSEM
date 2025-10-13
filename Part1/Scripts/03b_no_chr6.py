import pandas as pd
import numpy as np
import os

def load_and_filter(path, save_cleaned=True):
    df = pd.read_csv(path, sep='\t', engine='python')
    df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce')
    df = df[df['CHR'] != 6]
    if save_cleaned:
        base, ext = os.path.splitext(path)
        cleaned_path = f"{base}_no_chr6{ext}"
        df.to_csv(cleaned_path, index=False)
        print(f"Saved cleaned file to: {cleaned_path}")
    return df

if __name__ == '__main__':
    # path = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv"
    path = "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv"
    df = load_and_filter(path)
