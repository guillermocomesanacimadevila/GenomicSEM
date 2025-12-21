#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.pybin = "python3"

process qc {
    tag "${trait}_qc"

    input:
    tuple val(trait), val(qc_script), val(raw_in), val(qc_out), val(cases), val(controls)

    output:
    tuple val(trait), val(qc_out), val(cases), val(controls)

    script:
    """
    set -e
    WORKDIR=\$(pwd)
    LOGDIR=${workflow.launchDir}/logs/qc
    mkdir -p "\$LOGDIR"
    mkdir -p ${workflow.launchDir}/data/Main/${trait}/post-qc
    cd "${workflow.launchDir}"

    ${params.pybin} ${workflow.launchDir}/src/qc/${qc_script} \
      --in  ${workflow.launchDir}/data/Main/${trait}/${raw_in} \
      --out ${workflow.launchDir}/data/Main/${trait}/post-qc/${qc_out}

    cd "\$WORKDIR"
    touch ${trait}_qc.done
    cp ${trait}_qc.done "\$LOGDIR/"
    """
}

process neff {
    tag "${trait}_neff"

    input:
    tuple val(trait), val(qc_out), val(cases), val(controls)

    output:
    path "${trait}_neff.done"

    script:
    """
    set -e
    WORKDIR=\$(pwd)
    LOGDIR=${workflow.launchDir}/logs/qc
    mkdir -p "\$LOGDIR"
    cd "${workflow.launchDir}"

    ${params.pybin} ${workflow.launchDir}/src/qc/02_compute_neff.py \
      --in  ${workflow.launchDir}/data/Main/${trait}/post-qc/${qc_out} \
      --out ${workflow.launchDir}/data/Main/${trait}/post-qc/${trait}.ldsc_ready_neff.tsv \
      --cases ${cases} \
      --controls ${controls}

    cd "\$WORKDIR"
    touch ${trait}_neff.done
    cp ${trait}_neff.done "\$LOGDIR/"
    """
}

workflow {
    
    trait_info = Channel.of(
        tuple('AD',  '01a_exploratory_ad.py', 'Kunkle_etal_2019_IGAP_Summary_statistics_published.tsv', 'Kunkle_etal_2019_IGAP_Summary_statistics_published_ldsc_ready.tsv', 35274, 59163),
        tuple('SCZ', '01b_exploratory_sz.py', 'PGC3_SCZ_wave3.cleaned.tsv', 'PGC3_SCZ_wave3.cleaned_ldsc_ready.tsv', 67390, 94015),
        tuple('LON', '01c_exploratory_longevity.py','timmers2020_healthspan_lifespan_longevity.tsv', 'AGE_ldsc_ready.tsv', 354854, 354855)
    )

    qc_out = qc(trait_info)
    neff(qc_out)
}
