#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ad_bip {
  output:
    path "ad_bip_ldsc.done"
  script:
  """
  set -e
  cd "${workflow.launchDir}"
  Rscript scr/global_rg/ldsc_ad_bip.R
  mkdir -p Data/AD/post-ldsc Data/BIP/post-ldsc Data/logs
  mv AD.sumstats* Data/AD/post-ldsc/ 2>/dev/null || true
  mv BIP.sumstats* Data/BIP/post-ldsc/ 2>/dev/null || true
  mv *.log Data/logs/ 2>/dev/null || true
  cd -
  touch ad_bip_ldsc.done
  """
}

process scz_bip {
  output:
    path "scz_bip_ldsc.done"
  script:
  """
  set -e
  cd "${workflow.launchDir}"
  Rscript scr/global_rg/ldsc_scz_bip.R
  mkdir -p Data/SCZ/post-ldsc Data/BIP/post-ldsc Data/logs
  mv SCZ.sumstats* Data/SCZ/post-ldsc/ 2>/dev/null || true
  mv BIP.sumstats* Data/BIP/post-ldsc/ 2>/dev/null || true
  mv *.log Data/logs/ 2>/dev/null || true
  cd -
  touch scz_bip_ldsc.done
  """
}

workflow {
  ad_bip()
  scz_bip()
}
