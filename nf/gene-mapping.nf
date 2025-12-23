#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.ensembl = "/Users/c24102394/ensemble/mart_export.txt"
params.pybin = "python3"

process MapGenes {
    tag "${pheno1_prefix}_${pheno2_prefix}_locus${loc}"

    publishDir "${workflow.launchDir}/outputs/gene-mappings/locus_${loc}", mode: 'copy'

    input:
    tuple \
        val(pheno1_prefix), \
        val(pheno2_prefix), \
        val(pheno3_prefix), \
        val(loc), \
        path(pheno1_snps, stageAs: "pheno1_snps.txt"), \
        path(pheno1_eqtl, stageAs: "pheno1_eqtl.txt"), \
        path(pheno1_ci, stageAs: "pheno1_ci.txt"), \
        path(pheno2_snps, stageAs: "pheno2_snps.txt"), \
        path(pheno2_eqtl, stageAs: "pheno2_eqtl.txt"), \
        path(pheno2_ci, stageAs: "pheno2_ci.txt"), \
        path(pheno3_snps, stageAs: "pheno3_snps.txt"), \
        path(pheno3_eqtl, stageAs: "pheno3_eqtl.txt"), \
        path(pheno3_ci, stageAs: "pheno3_ci.txt"), \
        val(do_triple)

    output:
    path "*.tsv", emit: mappings
    path "${pheno1_prefix}_${pheno2_prefix}_${pheno3_prefix}_locus${loc}.done", emit: done

    script:
    """
    set -e
    mkdir -p ${workflow.launchDir}/logs/gene-maps

    ${params.pybin} ${workflow.launchDir}/src/map-genes/map.py \
        --pheno1_id ${pheno1_prefix} \
        --pheno2_id ${pheno2_prefix} \
        --pheno3_id ${pheno3_prefix} \
        --pheno1_snps pheno1_snps.txt \
        --pheno1_eqtl pheno1_eqtl.txt \
        --pheno1_ci pheno1_ci.txt \
        --pheno2_snps pheno2_snps.txt \
        --pheno2_eqtl pheno2_eqtl.txt \
        --pheno2_ci pheno2_ci.txt \
        --pheno3_snps pheno3_snps.txt \
        --pheno3_eqtl pheno3_eqtl.txt \
        --pheno3_ci pheno3_ci.txt \
        --ensembl_ref ${params.ensembl} \
        --out_dir . \
        \$( [[ "${do_triple}" == "true" ]] && echo "--do_triple" )

    touch ${pheno1_prefix}_${pheno2_prefix}_${pheno3_prefix}_locus${loc}.done
    cp .command.* ${workflow.launchDir}/logs/gene-maps/ || true
    """
}


// process FILTER_BY_SIGNIFICANC_AND_MULTIPLE_TESTING {
// 
// }

workflow {

    def fuma = file("${workflow.launchDir}/outputs/fuma_post")
    def pheno1 = 'AD'
    def pheno2 = 'SCZ'
    def pheno3 = 'LON'
    Channel.of(0, 1, 2)
        .map { loc ->
            tuple(
                pheno1, pheno2, pheno3, loc,
                fuma / "${pheno1}/locus_${loc}/snps.txt",
                fuma / "${pheno1}/locus_${loc}/eqtl.txt",
                fuma / "${pheno1}/locus_${loc}/ci.txt",
                fuma / "${pheno2}/locus_${loc}/snps.txt",
                fuma / "${pheno2}/locus_${loc}/eqtl.txt",
                fuma / "${pheno2}/locus_${loc}/ci.txt",
                fuma / "${pheno3}/locus_0/snps.txt",
                fuma / "${pheno3}/locus_0/eqtl.txt",
                fuma / "${pheno3}/locus_0/ci.txt",
                (loc == 2)
            )
        } \
        | MapGenes
}
