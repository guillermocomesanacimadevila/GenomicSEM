# Untangling Shared Genetic Components (AD -> SCZ)

> *Multi-layer, cross-trait discovery pipeline*

---

**Authors:**  
Guillermo Comesaña Cimadevila · Emily Simmonds · Dervis Salih · Nicholas Bray · Valentina Escott-Price 

---

[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![CI](https://img.shields.io/badge/CI-nextflow%20ready-blueviolet.svg)]()
[![Workflow](https://img.shields.io/badge/workflow-Nextflow-0DC09D.svg)]()
[![Containers](https://img.shields.io/badge/containers-Docker%20%7C%20Singularity-2496ED.svg)]()
[![Languages](https://img.shields.io/badge/languages-R%20%7C%20Python%20%7C%20Bash-276DC3.svg)]()
[![R](https://img.shields.io/badge/R-%E2%89%A54.2-276DC3.svg)]()
[![Python](https://img.shields.io/badge/python-%E2%89%A53.11-3776AB.svg)]()
[![Shell](https://img.shields.io/badge/shell-bash-4EAA25.svg)]()
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20WSL2-lightgrey.svg)]()
[![Reproducibility](https://img.shields.io/badge/reproducible-containers%20%2B%20lockfiles-informational.svg)]()
[![Docs](https://img.shields.io/badge/docs-auto_generated-green.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)]()

---


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

## Run GWAS-to-GWAS coloc!

```bash
Rscript coloc.R \
  AD_SCZ \
  AD \
  SCZ \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/loci_sumstats \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/coloc/AD_SCZ \
  cc cc \
  0.34 0.38
```
```bash
Rscript coloc.R \
  AD_LONG_SCZ_LONG \
  AD \
  LONG \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/loci_sumstats \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/coloc/AD_LONG \
  cc quant \
  0.34 0
```

```bash
Rscript coloc.R \
  AD_LONG_SCZ_LONG \
  SCZ \
  LONG \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/loci_sumstats \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/coloc/SCZ_LONG \
  cc quant \
  0.38 0
```


## Run MR! - You must have an embedded OpenGWAS token in your Renv

```bash
nextflow mr.nf --ld_token 
```

```bash
Rscript mr-pipeline.R \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv \
  /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/AGE/post-qc/timmers2020_healthspan_lifespan_longevity_neff.tsv \
  SCZ \
  LON \
  5e-8 \
  "NA"
```


## Run Mappings! 

```bash
python map.py \
  --pheno1_id AD \
  --pheno2_id SCZ \
  --pheno3_id LON \
  --pheno1_snps /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/AD/locus_2/snps.txt \
  --pheno1_eqtl /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/AD/locus_2/eqtl.txt \
  --pheno1_ci /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/AD/locus_2/ci.txt \
  --pheno2_snps /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/SCZ/locus_2/snps.txt \
  --pheno2_eqtl /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/SCZ/locus_2/eqtl.txt \
  --pheno2_ci /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/SCZ/locus_2/ci.txt \
  --pheno3_snps /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/LON/locus_0/snps.txt \
  --pheno3_eqtl /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/LON/locus_0/eqtl.txt \
  --pheno3_ci /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post/LON/locus_0/ci.txt \
  --ensembl_ref /Users/c24102394/ensemble/mart_export.txt \
  --out_dir /Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/gene-mappings/locus_2 \
  --do_triple
```


