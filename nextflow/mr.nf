#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.ad = '/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv'
params.scz = '/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv'
params.rbin = 'Rscript'
params.p_threshold = '5e-8'
params.ld_token = ''
params.conf_file = "${workflow.launchDir}/Data/MR/confounders.txt"

process RunMR4_AD_SCZ {
    publishDir "${workflow.launchDir}/outputs/MR/ad_scz", mode: 'copy'

    output:
        path "*"

    script:
    """
    set -euo pipefail
    ${params.rbin} ${workflow.launchDir}/scr/mr/run-mr.R \
        "${params.ad}" \
        "${params.scz}" \
        "AD" \
        "SCZ" \
        ${params.p_threshold} \
        "${params.ld_token}" \
        "${params.conf_file}"
    """
}

workflow {
    RunMR4_AD_SCZ()
}

// nextflow run mr.nf \
//   --ld_token 'LDLINK_TOKEN'
