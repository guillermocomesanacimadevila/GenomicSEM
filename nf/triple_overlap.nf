#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.pybin = 'python3'
params.thresh = 0.05

process TRIPLE_OVERLAP {
    tag "${out_name}"

    input:
    tuple val(conjfdr_dir), val(out_dir), val(out_name), val(phenos), val(thresh)

    output:
    path("overlap_snps.tsv")
    path("overlap_merged.tsv")

    script:
    """
    set -e
    ${params.pybin} ${workflow.launchDir}/src/pleio/make_triple_overlap.py \\
      --conjfdr-dir ${conjfdr_dir} \\
      --phenos ${phenos.join(' ')} \\
      --out-dir . \\
      --thresh ${thresh}

    mkdir -p ${out_dir}
    cp overlap_snps.tsv overlap_merged.tsv ${out_dir}/
    """
}

workflow {
    def conjfdr_dir = "${workflow.launchDir}/outputs/conjFDR"
    def out_name = "AD-SCZ-LON"
    def out_dir = "${workflow.launchDir}/outputs/conjFDR/${out_name}"
    def phenos = ["AD_LON", "SCZ_LON", "AD_SCZ"]

    Channel
        .of(tuple(conjfdr_dir, out_dir, out_name, phenos, params.thresh))
        | TRIPLE_OVERLAP
}
