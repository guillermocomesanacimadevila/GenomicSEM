#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ====================================== // 
// (BARE IN MIND -> ONLY 1,000 GENOMES)   //
// ====================================== //    

// baseline params 
// nextflow run sumher.nf --do_tagging false
// nextflow run sumher.nf --do_tagging true
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

    script:
    """
    set -euo pipefail
    mkdir -p ${workflow.launchDir}/logs/SumHer
    mkdir -p ${workflow.launchDir}/outputs/sumher
    mkdir -p ${workflow.launchDir}/outputs/sumher/tagging

    OUT="${workflow.launchDir}/outputs/sumher/tagging/eur_humdef"
    REF="${params.refpref}"

    for chr in {1..22}; do
        ${params.ldak} \
          --calc-tagging ${OUT}.chr${chr} \
          --bfile ${REF}.${chr} \
          --power ${params.power} \
          --chr ${chr}
    done

    ls ${OUT}.chr*.tagging > ${OUT}.taglist
    ${params.ldak} --join-tagging ${OUT} --taglist ${OUT}.taglist

    touch LDAK-tagging.done
    cp .command.* ${workflow.launchDir}/logs/SumHer/ || true

    """
}

// split h2 between AD and (SCZ & LON) -> different GWAS structure
process CALC_H2_AD {

    tag "${pheno_id}_SumHer_h2"

    input:
    tuple val(pheno_id), path(pheno_gwas)

    output:
    path "${pheno_id}_h2.*", emit: h2_all

    script:
    """
    set -euo pipefail
    mkdir -p ${workflow.launchDir}/outputs/sumher/${pheno_id}

    TAG="${workflow.launchDir}/outputs/sumher/tagging/eur_humdef.tagging"
    awk --version
    awk 'BEGIN{OFS="\t"} NR==1{print "Predictor","A1","A2","Z","N"} NR>1{print $1,$2,$3,$6/$7,$5}' \
        ${pheno_gwas} \
        > ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}.ldak.summaries
    
    # ldak
    ${params.ldak} \
        --sum-hers ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}_h2 \
        --summary ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}.ldak.summaries \
        --tagfile $TAG \
        --cutoff ${params.thresh} \
        --check-sums NO \
        > ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}_sumher.log 2>&1
    
    touch ${pheno_id}_SumHer_h2.done
    cp .command.* ${workflow.launchDir}/logs/SumHer/ || true
    ls -lh ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}_h2.hers
    head -5 ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}_h2.hers
    """
}

process CALC_H2_SCZ_LON {

    tag "${pheno_id}_SumHer_h2_B"

    input:
    tuple val(pheno_id), path(pheno_sumstats)

    output:
    path "${pheno_id}_h2.*", emit: h2_all_b

    script:
    """
    set -euo pipefail
    mkdir -p ${workflow.launchDir}/outputs/sumher/${pheno_id}
    
    TAG="${workflow.launchDir}/outputs/sumher/tagging/eur_humdef.tagging"
    awk 'BEGIN{OFS="\t"}
    NR==1{print "Predictor","A1","A2","Z","N"}
    NR>1{print $1,$2,$3,$6/$7,$5}' \
    ${pheno_sumstats} \
    > ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}.ldak.summaries

    ${params.ldak} \
        --sum-hers ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}_h2 \
        --summary ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}.ldak.summaries \
        --tagfile $TAG \
        --cutoff ${params.thresh} \
        --check-sums NO \
        > ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}_sumher.log 2>&1

    touch ${pheno_id}_SumHer_h2_B.done
    cp .command.* ${workflow.launchDir}/logs/SumHer/ || true
    ls -lh ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}_h2.hers
    head -5 ${workflow.launchDir}/outputs/sumher/${pheno_id}/${pheno_id}_h2.hers
    """
}

process CALC_RG {

    tag "${pheno1_id}_${pheno2_id}_SumHer_rg"

    input:
    tuple val(pheno1_prefix), val(pheno2_prefix)

    output:
    path "${pheno1_prefix}-${pheno2_prefix}.*", emit: rg_all

    script:
    """
    set -euo pipefail
    mkdir -p ${workflow.launchDir}/outputs/sumher/rg
    mkdir -p ${workflow.launchDir}/outputs/sumher/rg/${pheno1_prefix}-${pheno2_prefix}

    TAG="${workflow.launchDir}/outputs/sumher/tagging/eur_humdef.tagging"
    ${params.ldak}} \
        --sum-cors ${workflow.launchDir}/outputs/sumher/rg/${pheno1_prefix}-${pheno2_prefix}/${pheno1_prefix}-${pheno2_prefix} \
        --summary ${workflow.launchDir}/outputs/sumher/${pheno1_prefix}/${pheno1_prefix}.ldak.summaries \
        --summary2 ${workflow.launchDir}/outputs/sumher/${pheno2_prefix}/${pheno2_prefix}.ldak.summaries \
        --tagfile $TAG \
        --cutoff ${params.thresh} \
        --check-sums NO \
        > ${workflow.launchDir}/outputs/sumher/rg/${pheno1_prefix}-${pheno1_prefix}/${pheno1_prefix}-${pheno1_prefix}.log 2>&1
    
    touch ${pheno1_id}_${pheno2_id}_SumHer_rg.done
    cp .command.* ${workflow.launchDir}/logs/SumHer/ || true
    """
}

// process COMPUTE_PVALS {
// 
// }

workflow {

    
}
