#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.data_dir = 'Data'
params.rbin = 'Rscript'
params.pybin = 'python3'

process FormatGWAS {
    publishDir "${workflow.launchDir}/outputs/conjFDR", mode: 'copy'
    output:
        path "data/AD_SCZ.harmonised_AD_SCZ.tsv"

    script:
    """
    set -euo pipefail
    mkdir -p data
    ${params.pybin} ${workflow.launchDir}/scr/conjFDR_prep/conjFDR_prep.py \
      --ad ${workflow.launchDir}/Data/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv \
      --scz ${workflow.launchDir}/Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv \
      --out_prefix AD_SCZ
    mv ${workflow.launchDir}/scr/conjFDR/data/AD_SCZ.harmonised_AD_SCZ.tsv data/
    """
}

process RunConjFDR {
    publishDir "${workflow.launchDir}/outputs/conjFDR", mode: 'copy'
    input:
        path "data/AD_SCZ.harmonised_AD_SCZ.tsv"
    output:
        path "AD_SCZ_cfdr_results.tsv"
        path "AD_SCZ_shared_hits.tsv"

    script:
    """
    set -euo pipefail
    ${params.rbin} ${workflow.launchDir}/scr/conjFDR/conjFDR.R data/AD_SCZ.harmonised_AD_SCZ.tsv AD_SCZ
    mv ../../outputs/conjFDR/AD_SCZ_cfdr_results.tsv .
    mv ../../outputs/conjFDR/AD_SCZ_shared_hits.tsv .
    """
}

workflow {
    formatted = FormatGWAS()
    RunConjFDR(formatted)
}
