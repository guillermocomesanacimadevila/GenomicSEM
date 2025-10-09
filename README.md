# GenomicSEM

<img width="3200" height="2000" alt="two_trait_manhattan" src="https://github.com/user-attachments/assets/e9d67798-0c08-4aea-a4f0-dcc53c3fe220" />

## Manhattan plotting! 

¨¨¨bash
python manhattan.py \
  --trait1 "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv" \
  --trait2 "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv" \
  --out "/Users/guillermocomesanacimadevila/Desktop/two_trait_manhattan" \
  --genomewide 5e-8 \
  --suggestive 1e-5 \
  --width 12 \
  --height 8
  ¨¨¨
