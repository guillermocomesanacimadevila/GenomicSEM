#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.pybin = "python3"
params.rbin = "Rscript"
params.refdir = "/Users/c24102394/ref/conjFDR"

process prep_data {
    tag "${trait1}_${trait2}_prep"

    input:
    tuple val(trait1), val(trait2), path(trait1_file), path(trait2_file), val(out_prefix), val(out_dir)

    output:
    path("${trait1}_${trait2}_conjfdr_prep.done")

    script:
    """
    set -e
    mkdir -p ${out_dir}
    ${params.pybin} ${workflow.launchDir}/src/pleio/conjFDR_prep.py \\
        --trait1 ${trait1_file} \\
        --trait2 ${trait2_file} \\
        --prefix1 ${trait1} \\
        --prefix2 ${trait2} \\
        --out_prefix ${out_prefix} \\
        --out_dir ${out_dir}
    touch ${trait1}_${trait2}_conjfdr_prep.done
    mkdir -p ${workflow.launchDir}/logs/conjFDR
    cp ${trait1}_${trait2}_conjfdr_prep.done ${workflow.launchDir}/logs/conjFDR/
    """
}

process conjFDR {

    tag "${pheno_prefix}"

    input:
    tuple path(input_file),
          val(pheno_prefix),
          val(conjFDR_ref),
          val(out_dir)

    output:
    path "${pheno_prefix}_conjFDR.done", emit: conjFDR_done

    script:
    """
    set -e
    mkdir -p ${out_dir}
    mkdir -p ${workflow.launchDir}/logs/conjFDR

    ${params.rbin} ${workflow.launchDir}/src/pleio/conjFDR.R \\
        ${input_file} \\
        ${pheno_prefix} \\
        ${conjFDR_ref} \\
        ${out_dir}

    touch ${pheno_prefix}_conjFDR.done
    cp .command.* ${workflow.launchDir}/logs/conjFDR/ || true
    """
}

workflow {

    def base = "${workflow.launchDir}/data/Main"
    def cfdr_out_base = "${workflow.launchDir}/outputs/conjFDR"

    trait_info_prep = Channel.of(
        tuple(
            "AD",
            "SCZ",
            file("${base}/AD/post-qc/AD.ldsc_ready_neff.tsv"),
            file("${base}/SCZ/post-qc/SCZ.ldsc_ready_neff.tsv"),
            "AD_SCZ",
            "${base}/conjFDR/AD_SCZ"
        ),
        tuple(
            "AD",
            "LON",
            file("${base}/AD/post-qc/AD.ldsc_ready_neff.tsv"),
            file("${base}/LON/post-qc/LON.ldsc_ready_neff.tsv"),
            "AD_LON",
            "${base}/conjFDR/AD_LON"
        ),
        tuple(
            "SCZ",
            "LON",
            file("${base}/SCZ/post-qc/SCZ.ldsc_ready_neff.tsv"),
            file("${base}/LON/post-qc/LON.ldsc_ready_neff.tsv"),
            "SCZ_LON",
            "${base}/conjFDR/SCZ_LON"
        )
    )

    prep_data(trait_info_prep)

    cfdr_inputs = Channel.of(
        tuple(
            file("${base}/conjFDR/AD_SCZ/AD_SCZ.harmonised_AD_SCZ.tsv"),
            "AD_SCZ",
            params.refdir,
            "${cfdr_out_base}/AD_SCZ"
        ),
        tuple(
            file("${base}/conjFDR/AD_LON/AD_LON.harmonised_AD_LON.tsv"),
            "AD_LON",
            params.refdir,
            "${cfdr_out_base}/AD_LON"
        ),
        tuple(
            file("${base}/conjFDR/SCZ_LON/SCZ_LON.harmonised_SCZ_LON.tsv"),
            "SCZ_LON",
            params.refdir,
            "${cfdr_out_base}/SCZ_LON"
        )
    )

    conjFDR(cfdr_inputs)
}
