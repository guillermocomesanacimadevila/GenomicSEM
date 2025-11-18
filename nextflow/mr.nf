#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.ad = '/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv'
params.scz = '/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/Data/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv'
params.rbin = 'Rscript'
params.p_threshold = '5e-8'
params.ld_token = ''
params.conf_file = "${workflow.launchDir}/Data/MR/confounders_AD_SCZ.txt"

process CheckOpenGwas {
    env.OPENGWAS_JWT = System.getenv('OPENGWAS_JWT')

    output:
        path "check_opengwas_done"

    script:
    """
    set -euo pipefail

    if [ -z "\${OPENGWAS_JWT:-}" ]; then
      echo "ERROR: OPENGWAS_JWT not set. Export OPENGWAS_JWT in your shell before running Nextflow." >&2
      exit 1
    fi

    echo "OpenGWAS token is set." > check_opengwas_done
    """
}

process RunMR4_AD_SCZ {
    publishDir "${workflow.launchDir}/outputs/MR/ad_scz", mode: 'copy'
    output:
        path "*"
    script:
    """
    set -euo pipefail

    mkdir -p ${workflow.launchDir}/Data/MR

    cat << 'EOF' > ${params.conf_file}
smoking
bmi
body mass index
obesity
alcohol
education
schooling
socioeconomic
income
blood pressure
hypertension
stroke
cholesterol
hdl
ldl
triglyceride
diabetes
glucose
insulin
hba1c
inflammation
crp
immune
EOF

    ${params.rbin} ${workflow.launchDir}/scr/mr/mr-pipeline.R \
        "${params.ad}" \
        "${params.scz}" \
        "AD" \
        "SCZ" \
        ${params.p_threshold} \
        "${params.ld_token}" \
        "${params.conf_file}"
    """
}

process RunMR4_AD_LONG {

}

process RunMR4_SCZ_LONG {

}

process RunMR4_AD_SMOKING {

}

process RunMR4_SCZ_SMOKING {

}

workflow {
    CheckOpenGwas()
    RunMR4_AD_SCZ()
}
