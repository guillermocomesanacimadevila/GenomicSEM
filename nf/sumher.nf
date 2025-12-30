#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ====================================== // 
// (BARE IN MIND -> ONLY 1,000 GENOMES)   //
// ====================================== //    

// baseline params 
// nextflow run nf/sumher.nf --do_tagging false
// nextflow run nf/sumher.nf --do_tagging true
params.outdir = "${workflow.launchDir}/outputs/sumher"
params.pybin = "python3"
params.ldak = "/Users/c24102394/SumHer/LDAK/ldak6.1.mac" // HARD CODED -> change if you´re running nf script 
params.refpref = "/Users/c24102394/ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC" // HARD CODED -> change if you´re running nf script 
params.power = -0.25 // default param
params.thresh = 0.01 // default param
params.do_tagging = false // making sure that we dont re-create REFERENCE tagging files

process GEN_TAGGING_FILES {

    tag "LDAK-tagging"

    when:
    params.do_tagging

    output:
    path "eur_humdef.tagging"
    path "eur_humdef.taglist"

    publishDir "${workflow.launchDir}/outputs/sumher/tagging", mode: 'copy'

    script:
    """
    set -euo pipefail
    mkdir -p ${workflow.launchDir}/logs/SumHer

    OUT="eur_humdef"
    REF="${params.refpref}"

    for chr in {1..22}; do
        ${params.ldak} \\
          --calc-tagging ${OUT}.chr${chr} \\
          --bfile ${REF}.${chr} \\
          --power ${params.power} \\
          --chr ${chr}
    done

    ls ${OUT}.chr*.tagging > ${OUT}.taglist
    ${params.ldak} --join-tagging ${OUT} --taglist ${OUT}.taglist

    cp .command.* ${workflow.launchDir}/logs/SumHer/ || true
    """
}

// split h2 between AD and (SCZ & LON) -> different GWAS structure
process CALC_H2_AD {

    tag "${pheno_id}_SumHer_h2"

    input:
    tuple val(pheno_id), path(pheno_gwas)

    output:
    path "${pheno_id}.ldak.summaries", emit: sums
    path "${pheno_id}_h2.*", emit: h2
    path "${pheno_id}_sumher.log", emit: log

    publishDir "${workflow.launchDir}/outputs/sumher/${pheno_id}", mode: 'copy'

    script:
    """
    set -euo pipefail

    TAG="${workflow.launchDir}/outputs/sumher/tagging/eur_humdef.tagging"

    awk -v OFS='\\t' 'NR==1{print "Predictor","A1","A2","Z","N"; next} {print \$1,\$2,\$3,\$6/\$7,\$5}' \\
        ${pheno_gwas} > ${pheno_id}.ldak.summaries
    
    ${params.ldak} \\
        --sum-hers ${pheno_id}_h2 \\
        --summary ${pheno_id}.ldak.summaries \\
        --tagfile "\$TAG" \\
        --cutoff ${params.thresh} \\
        --check-sums NO \\
        > ${pheno_id}_sumher.log 2>&1

    cp .command.* ${workflow.launchDir}/logs/SumHer/ || true
    """
}

process CALC_H2_SCZ_LON {

    tag "${pheno_id}_SumHer_h2_B"

    input:
    tuple val(pheno_id), path(pheno_sumstats)

    output:
    path "${pheno_id}.ldak.summaries", emit: sums
    path "${pheno_id}_h2.*", emit: h2
    path "${pheno_id}_sumher.log", emit: log

    publishDir "${workflow.launchDir}/outputs/sumher/${pheno_id}", mode: 'copy'

    script:
    """
    set -euo pipefail

    TAG="${workflow.launchDir}/outputs/sumher/tagging/eur_humdef.tagging"

    awk -v OFS='\\t' 'NR==1{print "Predictor","A1","A2","Z","N"; next} {print \$1,\$2,\$3,\$6/\$7,\$5}' \\
        ${pheno_sumstats} > ${pheno_id}.ldak.summaries

    ${params.ldak} \\
        --sum-hers ${pheno_id}_h2 \\
        --summary ${pheno_id}.ldak.summaries \\
        --tagfile "\$TAG" \\
        --cutoff ${params.thresh} \\
        --check-sums NO \\
        > ${pheno_id}_sumher.log 2>&1

    cp .command.* ${workflow.launchDir}/logs/SumHer/ || true
    """
}

process CALC_RG {

    tag "${pheno1}-${pheno2}_SumHer_rg"

    input:
    tuple val(pheno1), val(pheno2)

    output:
    path "${pheno1}-${pheno2}.*", emit: rg_all
    path "res/*", emit: rg_res

    publishDir "${workflow.launchDir}/outputs/sumher/rg/${pheno1}-${pheno2}", mode: 'copy'

    script:
    """
    set -euo pipefail

    TAG="${workflow.launchDir}/outputs/sumher/tagging/eur_humdef.tagging"
    mkdir -p res

    ${params.ldak} \\
      --sum-cors "${pheno1}-${pheno2}" \\
      --summary "${workflow.launchDir}/outputs/sumher/${pheno1}/${pheno1}.ldak.summaries" \\
      --summary2 "${workflow.launchDir}/outputs/sumher/${pheno2}/${pheno2}.ldak.summaries" \\
      --tagfile "\$TAG" \\
      --cutoff ${params.thresh} \\
      --check-sums NO \\
      > "${pheno1}-${pheno2}.log" 2>&1

    ${params.pybin} ${workflow.launchDir}/src/sumher/calc_p.py \\
        --pheno1_prefix ${pheno1} \\
        --pheno2_prefix ${pheno2} \\
        --rg_column_name Value \\
        --se_column_name SE \\
        --results "${pheno1}-${pheno2}.cors" \\
        --out_dir res

    cp .command.* ${workflow.launchDir}/logs/SumHer/ || true
    """
}

// add pvalue calc per rg and oper trait pair save onto /results -> .tsv with       rg    P    SE

workflow {

    def base = "${workflow.launchDir}/data/Main"

    GEN_TAGGING_FILES()

    def ad_h2 = Channel.of(
        tuple(
            "AD",
            file("${base}/AD/post-qc/AD.ldsc_ready_neff.tsv")
        )
    )

    def scz_lon_h2 = Channel.of(
        tuple(
            "SCZ",
            file("${base}/SCZ/post-qc/SCZ.ldsc_ready_neff.tsv")
        ),
        tuple(
            "LON",
            file("${base}/LON/post-qc/LON.ldsc_ready_neff.tsv")
        )
    )

    CALC_H2_AD(ad_h2)
    CALC_H2_SCZ_LON(scz_lon_h2)

    def rg_pairs = Channel.of(
        tuple(
            "AD",
            "SCZ"
        ),
        tuple(
            "AD",
            "LON"
        ),
        tuple(
            "SCZ",
            "LON"
        )
    )

    CALC_RG(rg_pairs)
}
