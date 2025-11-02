#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.data_dir = 'Data'
params.logs_dir = "${params.data_dir}/logs"
params.rbin = 'Rscript'

process GenomicSEM {
  errorStrategy 'finish'

  output:
    path "genomic_sem.done"
    
  script:
  """
  set -euo pipefail
  taskdir=\$PWD
  cd "${workflow.launchDir}"

  if ! ${params.rbin} scr/genomicsem/models.R ; then
    echo "[GenomicSEM process] R failed. Check Data/logs/." >&2
  fi

  mkdir -p ${workflow.launchDir}/outputs
  mkdir -p ${workflow.launchDir}/outputs/GenomicSEM

  mv *.log ${params.logs_dir}/ 2>/dev/null || true
  mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true
  touch "\${taskdir}/genomicsem.done"
  """
}

workflow {
    GenomicSEM()
}