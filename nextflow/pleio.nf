#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.data_dir = 'Data'
params.rbin     = 'Rscript'
params.pybin    = 'python3'
params.ref_dir  = '/Users/c24102394/ref/ldsc/1000G_EUR_Phase3_plink'

process FormatGWAS_AD_SCZ {
    publishDir "${workflow.launchDir}/outputs/conjFDR", mode: 'copy'

    output:
        path "data/AD_SCZ.harmonised_AD_SCZ.tsv"

    script:
    """
    set -euo pipefail
    mkdir -p data
    ${params.pybin} ${workflow.launchDir}/scr/conjFDR/conjFDR_prep.py \
      --ad  ${workflow.launchDir}/${params.data_dir}/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv \
      --scz ${workflow.launchDir}/${params.data_dir}/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv \
      --out_prefix AD_SCZ
    mv ${workflow.launchDir}/scr/conjFDR/data/AD_SCZ.harmonised_AD_SCZ.tsv data/
    """
}

process FormatGWAS_AD_LONG {
    publishDir "${workflow.launchDir}/outputs/conjFDR", mode: 'copy'

    output:
        path "data/AD_LONG.harmonised_AD_LONG.tsv"

    script:
    """
    set -euo pipefail
    mkdir -p data
    ${params.pybin} ${workflow.launchDir}/scr/conjFDR/conjFDR_prep.py \
      --ad  ${workflow.launchDir}/${params.data_dir}/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv \
      --scz ${workflow.launchDir}/${params.data_dir}/AGE/post-qc/timmers2020_healthspan_lifespan_longevity_neff.tsv \
      --out_prefix AD_LONG
    mv ${workflow.launchDir}/scr/conjFDR/data/AD_LONG.harmonised_AD_LONG.tsv data/
    """
}

process FormatGWAS_SCZ_LONG {
    publishDir "${workflow.launchDir}/outputs/conjFDR", mode: 'copy'

    output:
        path "data/SCZ_LONG.harmonised_SCZ_LONG.tsv"

    script:
    """
    set -euo pipefail
    mkdir -p data
    ${params.pybin} ${workflow.launchDir}/scr/conjFDR/conjFDR_prep.py \
      --ad  ${workflow.launchDir}/${params.data_dir}/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv \
      --scz ${workflow.launchDir}/${params.data_dir}/AGE/post-qc/timmers2020_healthspan_lifespan_longevity_neff.tsv \
      --out_prefix SCZ_LONG
    mv ${workflow.launchDir}/scr/conjFDR/data/SCZ_LONG.harmonised_SCZ_LONG.tsv data/
    """
}

process RunConjFDR_AD_SCZ {
    publishDir "${workflow.launchDir}/outputs/conjFDR", mode: 'copy'

    input:
        path harmonised_ad_scz

    output:
        path "AD_SCZ_cfdr_results.tsv"
        path "AD_SCZ_shared_hits.tsv"

    script:
    """
    set -euo pipefail
    ${params.rbin} ${workflow.launchDir}/scr/conjFDR/conjFDR.R ${harmonised_ad_scz} AD_SCZ
    mv ${workflow.launchDir}/outputs/conjFDR/AD_SCZ_cfdr_results.tsv .
    mv ${workflow.launchDir}/outputs/conjFDR/AD_SCZ_shared_hits.tsv .
    """
}

process RunConjFDR_AD_LONG {
    publishDir "${workflow.launchDir}/outputs/conjFDR", mode: 'copy'

    input:
        path harmonised_ad_long

    output:
        path "AD_LONG_cfdr_results.tsv"
        path "AD_LONG_shared_hits.tsv"

    script:
    """
    set -euo pipefail
    ${params.rbin} ${workflow.launchDir}/scr/conjFDR/conjFDR.R ${harmonised_ad_long} AD_LONG
    mv ${workflow.launchDir}/outputs/conjFDR/AD_LONG_cfdr_results.tsv .
    mv ${workflow.launchDir}/outputs/conjFDR/AD_LONG_shared_hits.tsv .
    """
}

process RunConjFDR_SCZ_LONG {
    publishDir "${workflow.launchDir}/outputs/conjFDR", mode: 'copy'

    input:
        path harmonised_scz_long

    output:
        path "SCZ_LONG_cfdr_results.tsv"
        path "SCZ_LONG_shared_hits.tsv"

    script:
    """
    set -euo pipefail
    ${params.rbin} ${workflow.launchDir}/scr/conjFDR/conjFDR.R ${harmonised_scz_long} SCZ_LONG
    mv ${workflow.launchDir}/outputs/conjFDR/SCZ_LONG_cfdr_results.tsv .
    mv ${workflow.launchDir}/outputs/conjFDR/SCZ_LONG_shared_hits.tsv .
    """
}

process MapHits {
    publishDir "${workflow.launchDir}/outputs/conjFDR/mapped", mode: 'copy'

    input:
        path ad_scz_hits
        path ad_long_hits
        path scz_long_hits

    output:
        path "AD_SCZ_shared_hits_mapped.tsv"
        path "AD_LONG_shared_hits_mapped.tsv"
        path "SCZ_LONG_shared_hits_mapped.tsv"
        path "shared_snps_across_trait_pairs_full.csv"
        path "shared_snps_across_trait_pairs.csv"

    script:
    """
    set -euo pipefail
    ${params.pybin} ${workflow.launchDir}/scr/conjFDR/map_hits.py \
      --ref-dir ${params.ref_dir} \
      --out-dir . \
      ${ad_scz_hits} \
      ${ad_long_hits} \
      ${scz_long_hits}
    """
}

workflow {
    ad_scz_ch = FormatGWAS_AD_SCZ()
    ad_long_ch = FormatGWAS_AD_LONG()
    scz_long_ch = FormatGWAS_SCZ_LONG()

    ad_scz_conj_ch = RunConjFDR_AD_SCZ(ad_scz_ch)
    ad_long_conj_ch = RunConjFDR_AD_LONG(ad_long_ch)
    scz_long_conj_ch = RunConjFDR_SCZ_LONG(scz_long_ch)

    ad_scz_hits_ch = ad_scz_conj_ch.filter { it.name == "AD_SCZ_shared_hits.tsv" }
    ad_long_hits_ch = ad_long_conj_ch.filter { it.name == "AD_LONG_shared_hits.tsv" }
    scz_long_hits_ch = scz_long_conj_ch.filter { it.name == "SCZ_LONG_shared_hits.tsv" }

    MapHits(ad_scz_hits_ch, ad_long_hits_ch, scz_long_hits_ch)
}

// Rscript /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/scr/conjFDR/conjFDR.R data/AD_LONG.harmonised_AD_SCZ.tsv AD_LONG
// Rscript /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/scr/conjFDR/conjFDR.R data/SCZ_LONG.harmonised_AD_SCZ.tsv SCZ_LONG

// python map_hits.py \
//   --ref-dir /Users/c24102394/ref/ldsc/1000G_EUR_Phase3_plink \
//   --out-dir /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/conjFDR/mapped \
//   /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/conjFDR/AD_SCZ_shared_hits.tsv \
//   /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/conjFDR/AD_LONG_shared_hits.tsv \
//   /Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/conjFDR/SCZ_LONG_shared_hits.tsv
