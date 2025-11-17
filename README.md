# Shared Genetic Components (AD -> SCZ -> BIP)

<img width="2800" height="1800" alt="manhattan_balanced_colour_realdata" src="https://github.com/user-attachments/assets/82a652ad-b1b5-43ce-9464-5ca81539f4c5" />



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

```bash
python3 02_compute_neff.py \
  --in ../../Data/AD/post-qc/Kunkle_etal_2019_IGAP_Summary_statistics_published_ldsc_ready.tsv \
  --out ../../Data/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv \
  --cases 35274 \
  --controls 59163
```

```bash
python3 02_compute_neff.py \
  --in ../../Data/SCZ/post-qc/PGC3_SCZ_wave3.cleaned_ldsc_ready.tsv \
  --out ../../Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv \
  --cases 67390 \
  --controls 94015
```
## Run Experiments!

```bash
git clone https://github.com/guillermocomesanacimadevila/GenomicSEM.git
```

```bash
cd GenomicSEM
```

```bash
nextflow run nextflow/main.nf 
```

## Run MR!

```bash
nextflow mr.nf --ld_token 'LDLINK_TOKEN'
```

```bash
Rscript run-mr.R \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv \
  AD \
  SCZ \
  5e-8 \
  "LDLINK_TOKEN" \
  confounders_ad_scz.txt
```

