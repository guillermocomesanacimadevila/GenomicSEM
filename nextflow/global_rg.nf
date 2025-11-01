#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.data_dir = 'Data'
params.logs_dir = "${params.data_dir}/logs"
params.rbin = 'Rscript'

process ad_bip {
  errorStrategy 'finish'

  output:
    path "ad_bip_ldsc.done"

  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"
  if ! ${params.rbin} scr/global_rg/ldsc_ad_bip.R ; then
    echo "[ad_bip] R failed (probably at munge -> BIP.sumstats). Check Data/logs/." >&2
  fi

  mkdir -p ${params.data_dir}/AD/post-ldsc
  mkdir -p ${params.data_dir}/BIP/post-ldsc
  mkdir -p ${params.logs_dir}
  mkdir -p outputs

  # R writes into . AND into Data/, so sweep both
  mv AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/  2>/dev/null || true
  mv ${params.data_dir}/AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/  2>/dev/null || true

  mv BIP.sumstats*.gz ${params.data_dir}/BIP/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/BIP.sumstats*.gz ${params.data_dir}/BIP/post-ldsc/ 2>/dev/null || true

  # logs from both places
  mv *.log ${params.logs_dir}/ 2>/dev/null || true
  mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

  # tell NF we are done
  touch "\${taskdir}/ad_bip_ldsc.done"
  """
}

process scz_bip {
  errorStrategy 'finish'

  output:
    path "scz_bip_ldsc.done"

  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"

  if ! ${params.rbin} scr/global_rg/ldsc_scz_bip.R ; then
    echo "[scz_bip] R failed (probably at munge -> SCZ/BIP). Check Data/logs/." >&2
  fi

  mkdir -p ${params.data_dir}/SCZ/post-ldsc
  mkdir -p ${params.data_dir}/BIP/post-ldsc
  mkdir -p ${params.logs_dir}
  mkdir -p outputs

  mv SCZ.sumstats*.gz ${params.data_dir}/SCZ/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/SCZ.sumstats*.gz ${params.data_dir}/SCZ/post-ldsc/ 2>/dev/null || true

  mv BIP.sumstats*.gz ${params.data_dir}/BIP/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/BIP.sumstats*.gz ${params.data_dir}/BIP/post-ldsc/ 2>/dev/null || true

  mv *.log ${params.logs_dir}/ 2>/dev/null || true
  mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

  touch "\${taskdir}/scz_bip_ldsc.done"
  """
}

process ad_scz {
  errorStrategy 'finish'

  output:
    path "ad_scz_ldsc.done"

  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"

  rm -f AD.sumstats AD.sumstats.gz SCZ.sumstats SCZ.sumstats.gz \
        ${params.data_dir}/AD.sumstats ${params.data_dir}/AD.sumstats.gz \
        ${params.data_dir}/SCZ.sumstats ${params.data_dir}/SCZ.sumstats.gz 2>/dev/null || true

  if ! ${params.rbin} scr/global_rg/ldsc_ad_scz.R ; then
    echo "[ad_scz] R failed (probably at munge -> AD.sumstats). Check Data/logs/." >&2
  fi

  mkdir -p ${params.data_dir}/AD/post-ldsc
  mkdir -p ${params.data_dir}/SCZ/post-ldsc
  mkdir -p ${params.logs_dir}
  mkdir -p outputs

  mv AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/  2>/dev/null || true
  mv ${params.data_dir}/AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/  2>/dev/null || true

  mv SCZ.sumstats*.gz ${params.data_dir}/SCZ/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/SCZ.sumstats*.gz ${params.data_dir}/SCZ/post-ldsc/ 2>/dev/null || true

  mv *.log ${params.logs_dir}/ 2>/dev/null || true
  mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

  touch "\${taskdir}/ad_scz_ldsc.done"
  """
}

workflow {
  ad_bip()
  scz_bip()
  ad_scz()
}
