#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.rbin = "Rscript"

process LDSC {
    tag "${trait1}_${trait2}_ldsc"

    input:
    tuple val(trait1), val(trait2), val(pheno1_file), val(pheno2_file), val(cases1), val(controls1), val(cases2), val(controls2), val(pp1), val(pp2)

    output:
    path "${trait1}_${trait2}_ldsc.done", emit: ldsc_done

    script:
    """
    set -e
    mkdir -p ${workflow.launchDir}/outputs/ldsc
    mkdir -p ${workflow.launchDir}/data/Main/${trait1}/post-ldsc
    mkdir -p ${workflow.launchDir}/data/Main/${trait2}/post-ldsc

    ${params.rbin} ${workflow.launchDir}/src/ldsc/ldsc.R \
        ${workflow.launchDir}/data/Main/${trait1}/post-qc/${pheno1_file} \
        ${workflow.launchDir}/data/Main/${trait2}/post-qc/${pheno2_file} \
        ${trait1} \
        ${trait2} \
        ${cases1} \
        ${controls1} \
        ${cases2} \
        ${controls2} \
        ${pp1} \
        ${pp2} \
        ${workflow.launchDir}/outputs/ldsc

    touch ${trait1}_${trait2}_ldsc.done
    mkdir -p ${workflow.launchDir}/logs/ldsc
    mkdir -p ${workflow.launchDir}/logs/ldsc/${trait1}_${trait2}_ldsc
    cp ${trait1}_${trait2}_ldsc.done ${workflow.launchDir}/logs/ldsc/${trait1}_${trait2}_ldsc/
    """
}

workflow {
    trait_info = Channel.of(
        tuple(
            "AD",
            "SCZ",
            "AD.ldsc_ready_neff.tsv",
            "SCZ.ldsc_ready_neff.tsv",
            21982,
            41944,
            67390,
            94015,
            0.07,
            0.01
        ),
        tuple(
            "AD",
            "LON",
            "AD.ldsc_ready_neff.tsv",
            "LON.ldsc_ready_neff.tsv",
            21982,
            41944,
            354854,
            354855,
            0.07,
            0.10
        ),
        tuple(
            "SCZ",
            "LON",
            "SCZ.ldsc_ready_neff.tsv",
            "LON.ldsc_ready_neff.tsv",
            67390,
            94015,
            354854,
            354855,
            0.01,
            0.10
        )
    )

    LDSC(trait_info)
}
