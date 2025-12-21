#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.out_base = "${workflow.launchDir}/outputs/fuma_prep"
params.log_base = "${workflow.launchDir}/logs/fuma_prep"
params.clump_base = "${workflow.launchDir}/outputs/clumping"

process PREP_4_FUMA {
  // mkdir -p logs/fuma_prep
  // mkdir -p outputs/fuma_prep
  // mkdir -p outputs/fuma_prep/AD-SCZ
  // mkdir -p outputs/fuma_prep/AD-SCZ-LON
  // copy from clumping outputs each respective locus for each respetcie pheno
  // mkdir -p output/fuma_results
  // gzip the locus 0/1 file per trait
  tag "${set_id}"

  publishDir "${params.out_base}", mode: 'copy'

  input:
  tuple val(set_id), val(clump_dir), val(traits)

  output:
  path("${set_id}")

  script:
  def traits_str = traits.join(' ')
  """
  set -e
  mkdir -p ${set_id}/loci
  mkdir -p ${set_id}/meta
  echo "set_id\ttrait\tlocus_id\tfile_gz" > ${set_id}/meta/manifest.tsv
  for t in ${traits_str}; do
    mkdir -p ${set_id}/loci/\$t
    for f in ${clump_dir}/locus_*_\$t.tsv; do
      [ -f "\$f" ] || continue
      bn=\$(basename "\$f")
      locus_id=\$(echo "\$bn" | cut -d'_' -f2)
      gzip -c "\$f" > ${set_id}/loci/\$t/\${bn}.gz
      echo "${set_id}\t\$t\t\$locus_id\t${set_id}/loci/\$t/\${bn}.gz" >> ${set_id}/meta/manifest.tsv
    done
  done
  [ -f ${clump_dir}/lead_snps.tsv ] && cp ${clump_dir}/lead_snps.tsv ${set_id}/meta/
  [ -f ${clump_dir}/locus_coords.tsv ] && cp ${clump_dir}/locus_coords.tsv ${set_id}/meta/
  """
}

workflow {
  Channel.of(
    tuple("AD-SCZ", "${params.clump_base}/AD_SCZ", ["AD","SCZ"]),
    tuple("AD-SCZ-LON", "${params.clump_base}/AD-SCZ-LON", ["AD","SCZ","LON"])
  ) | PREP_4_FUMA | view
}