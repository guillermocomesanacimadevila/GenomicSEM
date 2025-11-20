#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ad_qc {
  tag 'ad_qc'
  output:
    path "ad_qc.done"
  script:
  """
  set -e
  WORKDIR=\$(pwd)
  cd "${workflow.launchDir}"
  python3 scr/qc/01a_exploratory_ad.py \
    --in Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.tsv \
    --out Data/AD/post-qc/Kunkle_etal_2019_IGAP_Summary_statistics_published_ldsc_ready.tsv
  echo "AD QC done"
  cd "\$WORKDIR"
  touch ad_qc.done
  """
}

process ad_neff {
  tag 'ad_neff'
  input:
    path ad_done
  script:
  """
  set -e
  cd "${workflow.launchDir}"
  python3 scr/qc/02_compute_neff.py \
    --in Data/AD/post-qc/Kunkle_etal_2019_IGAP_Summary_statistics_published_ldsc_ready.tsv \
    --out Data/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv \
    --cases 35274 \
    --controls 59163
  echo "AD neff done"
  """
}

process sz_qc {
  tag 'sz_qc'
  output:
    path "sz_qc.done"
  script:
  """
  set -e
  WORKDIR=\$(pwd)
  cd "${workflow.launchDir}"
  python3 scr/qc/01b_exploratory_sz.py \
    --in Data/SCZ/PGC3_SCZ_wave3.cleaned.tsv \
    --out Data/SCZ/post-qc/PGC3_SCZ_wave3.cleaned_ldsc_ready.tsv
  echo "SCZ QC done"
  cd "\$WORKDIR"
  touch sz_qc.done
  """
}

process scz_neff {
  tag 'scz_neff'
  input:
    path sz_done
  script:
  """
  set -e
  cd "${workflow.launchDir}"
  python3 scr/qc/02_compute_neff.py \
    --in Data/SCZ/post-qc/PGC3_SCZ_wave3.cleaned_ldsc_ready.tsv \
    --out Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv \
    --cases 67390 \
    --controls 94015
  echo "SCZ neff done"
  """
}

process bip_qc {
  tag 'bip_qc'
  output:
    path "bip_qc.done"
  script:
  """
  set -e
  WORKDIR=\$(pwd)
  cd "${workflow.launchDir}"
  python3 scr/qc/01d_exploratory_bipolar.py \
    --in Data/BIP/bip2024_eur_no23andMe \
    --out Data/BIP/post-qc/bip2024_eur_no23andMe_ldsc_ready.tsv
  echo "BIP QC done"
  cd "\$WORKDIR"
  touch bip_qc.done
  """
}

process bip_neff {
  tag 'bip_neff'
  input:
    path bip_done
  script:
  """
  set -e
  cd "${workflow.launchDir}"
  python3 scr/qc/02_compute_neff.py \
    --in Data/BIP/post-qc/bip2024_eur_no23andMe_ldsc_ready.tsv \
    --out Data/BIP/post-qc/bip2024_eur_no23andMe_ldsc_ready_neff.tsv \
    --cases 130064 \
    --controls 2301519
  echo "BIP neff done"
  """
}

process long_qc {
  tag 'long_qc'
  output:
    path "long_qc.done"
  script:
  """
  set -e
  WORKDIR=\$(pwd)
  cd "${workflow.launchDir}"
  python3 scr/qc/01c_exploratory_longevity.py \
    --in Data/AGE/timmers2020_healthspan_lifespan_longevity.tsv \
    --out Data/AGE/post-qc/AGE_ldsc_ready.tsv
  echo "AGE QC done"
  cd "\$WORKDIR"
  touch long_qc.done
  """
}

process long_neff {
  tag 'long_neff'
  input:
    path long_done
  script:
  """
  set -e
  cd "${workflow.launchDir}"
  python3 scr/qc/02_compute_neff.py \
    --in Data/AGE/post-qc/AGE_ldsc_ready.tsv \
    --out Data/AGE/post-qc/timmers2020_healthspan_lifespan_longevity_neff.tsv \
    --cases 354854 \
    --controls 354855
  echo "LONG neff done"
  """
}

process smoking_qc {
  tag 'smoking_qc'
  output:
    path "smoking_qc.done"
  script:
  """
  set -e
  WORKDIR=\$(pwd)
  cd "${workflow.launchDir}"
  python3 scr/qc/01e_exploratory_smoking.py \
    --in Data/SMOKING/ieu-a-964.vcf \
    --out Data/SMOKING/post-qc/SMK_ldsc_ready.tsv
  echo "SMOKING QC done"
  cd "\$WORKDIR"
  touch smoking_qc.done
  """
}

process smoking_neff { // Neff = N -> 47,961
  tag 'smoking_neff'
  input:
    path smoking_neff
  script:
  """
  set -e
  cd "${workflow.launchDir}"
  python3 scr/qc/02_compute_neff.py \
    --in Data/SMOKING/post-qc/SMK_ldsc_ready.tsv \
    --out Data/SMOKING/post-qc/SMK_neff.tsv \
    --cases 23981 \
    --controls 23980
  echo "SMK neff done"
  """
}

workflow {
  ad_done = ad_qc()
  sz_done = sz_qc()
  bip_done = bip_qc()
  long_done = long_qc()
  smoking_done = smoking_qc()

  ad_neff(ad_done)
  scz_neff(sz_done)
  bip_neff(bip_done)
  long_neff(long_done)
  smoking_neff(smoking_done)
}
