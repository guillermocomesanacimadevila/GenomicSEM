#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.data_dir = 'Data'
params.logs_dir = "${params.data_dir}/logs"
params.pybin    = 'python3'

process prep_ad {
    errorStrategy 'finish'

    output:
    path "ad_lava_prep.done"

    script:
    """
    set -euo pipefail

    taskdir=\$PWD
    cd "${workflow.launchDir}"

    mkdir -p ${params.data_dir}/AD/post-lava
    mkdir -p ${params.data_dir}/AD/post-ldsc
    mkdir -p ${params.logs_dir}
    mkdir -p outputs

    ${params.pybin} scr/local_rg/prep_data.py \\
        -i ${params.data_dir}/AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv \\
        -o ${params.data_dir}/AD/post-lava/AD_lava_ready.tsv

    mv *.log ${params.logs_dir}/ 2>/dev/null || true
    mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

    touch "\${taskdir}/ad_lava_prep.done"
    """
}

process prep_scz {
    errorStrategy 'finish'

    output:
    path "scz_lava_prep.done"

    script:
    """
    set -euo pipefail

    taskdir=\$PWD
    cd "${workflow.launchDir}"

    mkdir -p ${params.data_dir}/SCZ/post-lava
    mkdir -p ${params.data_dir}/SCZ/post-ldsc
    mkdir -p ${params.logs_dir}
    mkdir -p outputs

    ${params.pybin} scr/local_rg/prep_data.py \\
        -i ${params.data_dir}/SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv \\
        -o ${params.data_dir}/SCZ/post-lava/SCZ_lava_ready.tsv

    mv *.log ${params.logs_dir}/ 2>/dev/null || true
    mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

    touch "\${taskdir}/scz_lava_prep.done"
    """
}

process prep_bip {
    errorStrategy 'finish'

    output:
    path "bip_lava_prep.done"

    script:
    """
    set -euo pipefail

    taskdir=\$PWD
    cd "${workflow.launchDir}"

    mkdir -p ${params.data_dir}/BIP/post-lava
    mkdir -p ${params.data_dir}/BIP/post-ldsc
    mkdir -p ${params.logs_dir}
    mkdir -p outputs

    ${params.pybin} scr/local_rg/prep_data.py \\
        -i ${params.data_dir}/BIP/post-qc/bip2024_eur_no23andMe_ldsc_ready_neff.tsv \\
        -o ${params.data_dir}/BIP/post-lava/BIP_lava_ready.tsv

    mv *.log ${params.logs_dir}/ 2>/dev/null || true
    mv ${params.data_dir}/*.log ${params.logs_dir}/ 2>/dev/null || true

    touch "\${taskdir}/bip_lava_prep.done"
    """
}

process make_input_file {
    errorStrategy 'finish'

    input:
    path ad_done
    path scz_done
    path bip_done

    output:
    path "InputFiles"
    path "input_files.done"

    script:
    """
    set -euo pipefail

    mkdir -p InputFiles

    # AD–SCZ
    cat > InputFiles/input.info_ad_scz.txt <<EOF
phenotype	cases	controls	prevalence	filename
AD	35274	59163	0.07	${workflow.launchDir}/Data/AD/post-lava/AD_lava_ready.tsv
SCZ	67390	94015	0.01	${workflow.launchDir}/Data/SCZ/post-lava/SCZ_lava_ready.tsv
EOF

    # SCZ–BIP
    cat > InputFiles/input.info_scz_bip.txt <<EOF
phenotype	cases	controls	prevalence	filename
SCZ	67390	94015	0.01	${workflow.launchDir}/Data/SCZ/post-lava/SCZ_lava_ready.tsv
BIP	130064	2301519	0.01	${workflow.launchDir}/Data/BIP/post-lava/BIP_lava_ready.tsv
EOF

    # AD–BIP
    cat > InputFiles/input.info_ad_bip.txt <<EOF
phenotype	cases	controls	prevalence	filename
AD	35274	59163	0.07	${workflow.launchDir}/Data/AD/post-lava/AD_lava_ready.tsv
BIP	130064	2301519	0.01	${workflow.launchDir}/Data/BIP/post-lava/BIP_lava_ready.tsv
EOF

    touch input_files.done
    mkdir -p ${workflow.launchDir}/Data/logs
    cp input_files.done ${workflow.launchDir}/Data/logs/
    """
}


workflow {
    ch_ad  = prep_ad()
    ch_scz = prep_scz()
    ch_bip = prep_bip()

    make_input_file(ch_ad, ch_scz, ch_bip)
}
