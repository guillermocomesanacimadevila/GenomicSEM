#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.ad = '/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv'
params.scz = '/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv'
params.age = '/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/AGE/post-qc/timmers2020_healthspan_lifespan_longevity_neff.tsv'
params.prefix = 'AD_SCZ_AGE'
params.outdir_asset = "${workflow.launchDir}/outputs/ASSET/ad_scz_age"
params.outdir_summary = "${workflow.launchDir}/outputs/summary"
params.pybin = 'python3'
params.rbin = 'Rscript'
params.scripts = "${workflow.launchDir}/scr/pleio"

process PrepData3 {
  publishDir "${params.outdir_asset}", mode: "copy"
  output:
  path "${params.prefix}.ASSET_input.tsv"
  script:
  """
  set -euo pipefail
  mkdir -p data
  mkdir -p ${params.outdir_asset}
  ${params.pybin} ${params.scripts}/asset_prep_data.py \
    --t1 ${params.ad} --t1-prefix AD \
    --t2 ${params.scz} --t2-prefix SCZ \
    --t3 ${params.age} --t3-prefix AGE \
    --out ${params.prefix}.ASSET_input.tsv
  """
}

process RunAsset {
  publishDir "${params.outdir_asset}", mode: "copy"
  input:
  path asset_input
  output:
  path "${params.prefix}_ASSET_results.tsv"
  path "${params.prefix}_ASSET_pleiotropic_hits.tsv"
  script:
  """
  set -euo pipefail
  mkdir -p data
  mkdir -p ${params.outdir_asset}
  cp ${asset_input} data/${params.prefix}.ASSET_input.tsv
  cd ${params.scripts}
  ${params.rbin} other_Asset.R data/${params.prefix}.ASSET_input.tsv ${params.prefix}
  cp ${params.outdir_asset}/${params.prefix}_ASSET_results.tsv ${PWD}/..
  cp ${params.outdir_asset}/${params.prefix}_ASSET_pleiotropic_hits.tsv ${PWD}/..
  """
}

process LogFindings {
  publishDir "${params.outdir_summary}", mode: "copy"
  input:
  path results_tsv
  path hits_tsv
  output:
  path "${params.prefix}_summary.txt"
  script:
  """
  set -euo pipefail
  mkdir -p ${params.outdir_summary}
  ${params.pybin} ${params.scripts}/print_summary.py \
    --results ${results_tsv} \
    --hits ${hits_tsv} \
    --out ${params.prefix}_summary.txt
  """
}

workflow {
  ch_in = PrepData3()
  ch_out = RunAsset(ch_in)
  LogFindings(ch_out)
}
