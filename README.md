# GenomicSEM

<img width="3200" height="2000" alt="three_manhattan_rw" src="https://github.com/user-attachments/assets/465249c9-db1b-487b-b41c-d8b4d78cb03a" />



## Manhattan plotting! 

```bash
python 03_manhattan.py \
  --trait1 "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv" \
  --trait2 "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv" \
  --out "/Users/guillermocomesanacimadevila/Desktop/two_trait_manhattan" \
  --genomewide 5e-8 \
  --suggestive 1e-5 \
  --width 12 \
  --height 8
```

```bash
Rscript ../Scripts/plots.R \
  SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv \
  AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv \
  manhattan_full.pdf
```
