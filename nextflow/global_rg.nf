#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.data_dir = 'Data'
params.logs_dir = "${params.data_dir}/logs"
params.rbin     = 'Rscript'

process ad_scz {
  errorStrategy 'finish'
  output:
    path "ad_scz_ldsc.done"
  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"

  rm -f AD.sumstats* SCZ.sumstats* \\
        ${params.data_dir}/AD.sumstats* ${params.data_dir}/SCZ.sumstats* 2>/dev/null || true

  if ! ${params.rbin} scr/global_rg/ldsc_ad_scz.R ; then
    echo "[ad_scz] R failed" >&2
    exit 1
  fi

  mkdir -p ${params.data_dir}/AD/post-ldsc
  mkdir -p ${params.data_dir}/SCZ/post-ldsc
  mkdir -p ${params.logs_dir}
  mkdir -p outputs

  mv AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/ 2>/dev/null || true

  mv SCZ.sumstats*.gz ${params.data_dir}/SCZ/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/SCZ.sumstats*.gz ${params.data_dir}/SCZ/post-ldsc/ 2>/dev/null || true

  mv *.log ${params.logs_dir}/ 2>/dev/null || true
  mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

  touch "\${taskdir}/ad_scz_ldsc.done"
  """
}

process ad_bip {
  errorStrategy 'finish'
  input:
    path ad_scz_done
  output:
    path "ad_bip_ldsc.done"
  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"

  rm -f AD.sumstats* BIP.sumstats* \\
        ${params.data_dir}/AD.sumstats* ${params.data_dir}/BIP.sumstats* 2>/dev/null || true

  if ! ${params.rbin} scr/global_rg/ldsc_ad_bip.R ; then
    echo "[ad_bip] R failed" >&2
    exit 1
  fi

  mkdir -p ${params.data_dir}/AD/post-ldsc
  mkdir -p ${params.data_dir}/BIP/post-ldsc
  mkdir -p ${params.logs_dir}
  mkdir -p outputs

  mv AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/ 2>/dev/null || true

  mv BIP.sumstats*.gz ${params.data_dir}/BIP/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/BIP.sumstats*.gz ${params.data_dir}/BIP/post-ldsc/ 2>/dev/null || true

  mv *.log ${params.logs_dir}/ 2>/dev/null || true
  mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

  touch "\${taskdir}/ad_bip_ldsc.done"
  """
}

process scz_bip {
  errorStrategy 'finish'
  input:
    path ad_bip_done
  output:
    path "scz_bip_ldsc.done"
  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"

  rm -f SCZ.sumstats* BIP.sumstats* \\
        ${params.data_dir}/SCZ.sumstats* ${params.data_dir}/BIP.sumstats* 2>/dev/null || true

  if ! ${params.rbin} scr/global_rg/ldsc_scz_bip.R ; then
    echo "[scz_bip] R failed" >&2
    exit 1
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

process ad_long {
  errorStrategy 'finish'
  input:
    path scz_bip_done
  output:
    path "ad_long_ldsc.done"
  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"

  rm -f AD.sumstats* AGE.sumstats* \\
        ${params.data_dir}/AD.sumstats* ${params.data_dir}/AGE.sumstats* 2>/dev/null || true

  if ! ${params.rbin} scr/global_rg/ldsc_ad_long.R ; then
    echo "[ad_long] R failed" >&2
    exit 1
  fi

  mkdir -p ${params.data_dir}/AD/post-ldsc
  mkdir -p ${params.data_dir}/AGE/post-ldsc
  mkdir -p ${params.logs_dir}
  mkdir -p outputs

  mv AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/AD.sumstats*.gz ${params.data_dir}/AD/post-ldsc/ 2>/dev/null || true

  mv AGE.sumstats*.gz ${params.data_dir}/AGE/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/AGE.sumstats*.gz ${params.data_dir}/AGE/post-ldsc/ 2>/dev/null || true

  mv *.log ${params.logs_dir}/ 2>/dev/null || true
  mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

  touch "\${taskdir}/ad_long_ldsc.done"
  """
}

process scz_long {
  errorStrategy 'finish'
  input:
    path ad_long_done
  output:
    path "scz_long_ldsc.done"
  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"

  rm -f SCZ.sumstats* AGE.sumstats* \\
        ${params.data_dir}/SCZ.sumstats* ${params.data_dir}/AGE.sumstats* 2>/dev/null || true

  if ! ${params.rbin} scr/global_rg/ldsc_scz_long.R ; then
    echo "[scz_long] R failed" >&2
    exit 1
  fi

  mkdir -p ${params.data_dir}/SCZ/post-ldsc
  mkdir -p ${params.data_dir}/AGE/post-ldsc
  mkdir -p ${params.logs_dir}
  mkdir -p outputs

  mv SCZ.sumstats*.gz ${params.data_dir}/SCZ/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/SCZ.sumstats*.gz ${params.data_dir}/SCZ/post-ldsc/ 2>/dev/null || true

  mv AGE.sumstats*.gz ${params.data_dir}/AGE/post-ldsc/ 2>/dev/null || true
  mv ${params.data_dir}/AGE.sumstats*.gz ${params.data_dir}/AGE/post-ldsc/ 2>/dev/null || true

  mv *.log ${params.logs_dir}/ 2>/dev/null || true
  mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

  touch "\${taskdir}/scz_long_ldsc.done"
  """
}

workflow {
  ad_scz_done  = ad_scz()
  ad_bip_done  = ad_bip(ad_scz_done)
  scz_bip_done = scz_bip(ad_bip_done)
  ad_long_done = ad_long(scz_bip_done)
  scz_long(ad_long_done)
}
