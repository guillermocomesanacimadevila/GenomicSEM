#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// NEED to change this to wherever your files are 
params.g100eur_prefix = "/Users/c24102394/MAGMA/ref/g1000_eur"
params.g100eur_bim = "/Users/c24102394/MAGMA/ref/g1000_eur.bim"
params.grch37_gene_loc = "/Users/c24102394/MAGMA/ref/NCBI37.3.gene.loc"
params.magma_bin = "/Users/c24102394/MAGMA/magma"

process posMAGMA {
  tag "${pheno_id}_magma"

  input:
  tuple val(pheno_id), path(pheno_file)

  output:
  path("${pheno_id}_magma.genes.out"), emit: genes_out
  path("${pheno_id}_magma.genes.mapped.tsv"), emit: mapped
  path("${pheno_id}_for_magma.txt"), emit: snp_pvals
  path("${pheno_id}_01magma.done"), emit: done

  publishDir "${workflow.launchDir}/outputs/magma/${pheno_id}/MAGMA", mode: 'copy'

  script:
  """
  set -euo pipefail
  chmod +x "${workflow.launchDir}/src/magma/01magma.sh"
  MAGMA_BIN="${params.magma_bin}" "${workflow.launchDir}/src/magma/01magma.sh" \
    . \
    "${params.g100eur_prefix}" \
    "${params.g100eur_bim}" \
    "${params.grch37_gene_loc}" \
    "${pheno_id}" \
    "${pheno_file}"

  touch "${pheno_id}_01magma.done"
  mkdir -p "${workflow.launchDir}/logs/magma/${pheno_id}/MAGMA"
  cp .command.* "${workflow.launchDir}/logs/magma/${pheno_id}/MAGMA/" 2>/dev/null || true
  """
}

workflow {
  def base = "${workflow.launchDir}/data/Main"

  magma_info = Channel.of(
    tuple("AD",  file("${base}/AD/post-ldsc/AD.sumstats.gz")),
    tuple("SCZ", file("${base}/SCZ/post-ldsc/SCZ.sumstats.gz")),
    tuple("LON", file("${base}/LON/post-ldsc/LON.sumstats.gz"))
  )

  posMAGMA(magma_info)
}

// process eMAGMA {
// 
// }

// process H_MAGMA {
// 
// }
