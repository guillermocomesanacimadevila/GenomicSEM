#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.rbin = "Rscript"
params.coloc_r = "${workflow.launchDir}/src/gwas-coloc/coloc.R"
params.out_base = "${workflow.launchDir}/outputs/coloc"

process PAIRWISE_COLOC {
  tag "${prefix}_${trait1}_${trait2}"

  input:
  tuple val(prefix),
        val(trait1),
        val(trait2),
        val(loci_dir),
        val(out_dir),
        val(type1),
        val(type2),
        val(s1),
        val(s2)

  output:
  path("${prefix}_${trait1}_${trait2}_coloc.tsv")

  script:
  """
  set -e
  mkdir -p ${out_dir}
  mkdir -p ${workflow.launchDir}/logs/coloc/pairwise
  ${params.rbin} ${params.coloc_r} \\
    ${prefix} \\
    ${trait1} \\
    ${trait2} \\
    ${loci_dir} \\
    ${out_dir} \\
    ${type1} \\
    ${type2} \\
    ${s1} \\
    ${s2}
  cp ${out_dir}/${prefix}_${trait1}_${trait2}_coloc.tsv ./
  cp .command.* ${workflow.launchDir}/logs/coloc/pairwise/ || true
  """
}

process TRIPLE_OVERLAP_COLOC {
  tag "${prefix}_TRIPLE"

  input:
  tuple val(prefix),
        val(trait1),
        val(trait2),
        val(trait3),
        val(loci_dir),
        val(out_dir),
        val(s1),
        val(s2)

  output:
  path("${prefix}_${trait1}_${trait3}_coloc.tsv")
  path("${prefix}_${trait2}_${trait3}_coloc.tsv")

  script:
  """
  set -e
  mkdir -p ${out_dir}
  mkdir -p ${workflow.launchDir}/logs/coloc/triple
  ${params.rbin} ${params.coloc_r} \\
    ${prefix} \\
    ${trait1} \\
    ${trait3} \\
    ${loci_dir} \\
    ${out_dir} \\
    cc \\
    quant \\
    ${s1}

  ${params.rbin} ${params.coloc_r} \\
    ${prefix} \\
    ${trait2} \\
    ${trait3} \\
    ${loci_dir} \\
    ${out_dir} \\
    cc \\
    quant \\
    ${s2}
  cp ${out_dir}/${prefix}_${trait1}_${trait3}_coloc.tsv ./
  cp ${out_dir}/${prefix}_${trait2}_${trait3}_coloc.tsv ./
  cp .command.* ${workflow.launchDir}/logs/coloc/triple/ || true
  """
}

workflow {

  def clump_pair_dir = "${workflow.launchDir}/outputs/clumping/AD_SCZ"
  def clump_triple_dir = "${workflow.launchDir}/outputs/clumping/AD-SCZ-LON"
  def coloc_base = params.out_base
  def ad_cases = 35274
  def ad_ctrls = 59163
  def scz_cases = 67390
  def scz_ctrls = 94015
  def s_ad = ad_cases / (ad_cases + ad_ctrls)
  def s_scz = scz_cases / (scz_cases + scz_ctrls)

  Channel.of(
    tuple(
      "AD_SCZ",
      "AD",
      "SCZ",
      clump_pair_dir,
      "${coloc_base}/AD_SCZ",
      "cc",
      "cc",
      s_ad,
      s_scz
    )
  ) | PAIRWISE_COLOC

  Channel.of(
    tuple(
      "AD-SCZ-LON",
      "AD",
      "SCZ",
      "LON",
      clump_triple_dir,
      "${coloc_base}/AD-SCZ-LON",
      s_ad,
      s_scz
    )
  ) | TRIPLE_OVERLAP_COLOC
}
