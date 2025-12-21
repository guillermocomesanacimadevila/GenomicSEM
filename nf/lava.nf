#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.pybin = "python3"
params.rbin = "Rscript"

process data_prep {
    tag "${trait}_prep4_lava"

    input:
    tuple val(trait), val(input_pheno)

    output:
    path("${trait}_lava_prep.done")

    script:
    """
    set -e
    mkdir -p ${workflow.launchDir}/data/Main/${trait}/lava-ready
    ${params.pybin} ${workflow.launchDir}/src/lava/prep_data.py \\
        --input ${workflow.launchDir}/data/Main/${trait}/post-qc/${input_pheno} \\
        --output ${workflow.launchDir}/data/Main/${trait}/lava-ready/${trait}.tsv
    touch ${trait}_lava_prep.done
    mkdir -p ${workflow.launchDir}/logs/lava
    cp ${trait}_lava_prep.done ${workflow.launchDir}/logs/lava/
    """
}

process infoTXT {

    output:
    path "input.info_ad_scz_age.txt"

    script:
    """
    set -e
    mkdir -p ${workflow.launchDir}/InputFiles
    cat << 'EOF' > input.info_ad_scz_age.txt
phenotype	cases	controls	prevalence	filename
AD	21982	41944	0.07	data/Main/AD/lava-ready/AD.tsv
LON	354854	354855	0.1	data/Main/LON/lava-ready/LON.tsv
SCZ	67390	94015	0.01	data/Main/SCZ/lava-ready/SCZ.tsv
EOF
    cp input.info_ad_scz_age.txt ${workflow.launchDir}/InputFiles/
    """
}

process LAVA {

    input:
    tuple  val(lava_ref),
           val(loci_file),
           val(info_tsv),
           val(ov12_file),
           val(ov13_file),
           val(ov23_file),
           val(out_dir),
           val(ph1_name), val(ph1_cases), val(ph1_ctrls), val(ph1_prev), path(ph1_file),
           val(ph2_name), val(ph2_cases), val(ph2_ctrls), val(ph2_prev), path(ph2_file),
           val(ph3_name), val(ph3_cases), val(ph3_ctrls), val(ph3_prev), path(ph3_file)

    output:
    path "lava.done", emit: lava_done

    script:
    """
    set -e
    mkdir -p ${out_dir}
    mkdir -p ${workflow.launchDir}/logs/lava

    ${params.rbin} ${workflow.launchDir}/src/lava/lava.R \\
        ${lava_ref} \\
        ${loci_file} \\
        ${info_tsv} \\
        ${ov12_file} \\
        ${ov13_file} \\
        ${ov23_file} \\
        ${out_dir} \\
        ${ph1_name} ${ph1_cases} ${ph1_ctrls} ${ph1_prev} ${ph1_file} \\
        ${ph2_name} ${ph2_cases} ${ph2_ctrls} ${ph2_prev} ${ph2_file} \\
        ${ph3_name} ${ph3_cases} ${ph3_ctrls} ${ph3_prev} ${ph3_file}

    touch lava.done
    cp .command.* ${workflow.launchDir}/logs/lava/ || true
    """
}

workflow {

    trait_info_prep = Channel.of(
        tuple("AD",  "AD.ldsc_ready_neff.tsv"),
        tuple("SCZ", "SCZ.ldsc_ready_neff.tsv"),
        tuple("LON", "LON.ldsc_ready_neff.tsv")
    )

    data_prep(trait_info_prep)
    infoTXT()

    lava_args = Channel.of(
        tuple(
            "/Users/c24102394/ref/lava/lava_ref/lava-ukb-v1.1",
            "/Users/c24102394/ref/lava/hdll_blocks.coords.loci",
            "${workflow.launchDir}/InputFiles/input.info_ad_scz_age.txt",
            "${workflow.launchDir}/outputs/ldsc/AD_LON/overlap_corr_for_LAVA_AD_LON.csv",
            "${workflow.launchDir}/outputs/ldsc/AD_SCZ/overlap_corr_for_LAVA_AD_SCZ.csv",
            "${workflow.launchDir}/outputs/ldsc/SCZ_LON/overlap_corr_for_LAVA_SCZ_LON.csv",
            "${workflow.launchDir}/outputs/LAVA/AD_LON_SCZ",
            "AD", 21982, 41944, 0.07,
            file("${workflow.launchDir}/data/Main/AD/lava-ready/AD.tsv"),
            "LON", 354854, 354855, 0.10,
            file("${workflow.launchDir}/data/Main/LON/lava-ready/LON.tsv"),
            "SCZ", 67390, 94015, 0.01,
            file("${workflow.launchDir}/data/Main/SCZ/lava-ready/SCZ.tsv")
        )
    )

    LAVA(lava_args)
}
