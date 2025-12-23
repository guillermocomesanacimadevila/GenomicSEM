#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.nfbin = "nextflow run"

process RUN_QC { script: "set -euo pipefail; ${params.nfbin} ${workflow.launchDir}/nf/qc.nf" }
process RUN_LDSC { script: "set -euo pipefail; ${params.nfbin} ${workflow.launchDir}/nf/ldsc.nf" }
process RUN_LAVA { script: "set -euo pipefail; ${params.nfbin} ${workflow.launchDir}/nf/lava.nf" }
process RUN_CONJFDR { script: "set -euo pipefail; ${params.nfbin} ${workflow.launchDir}/nf/pleio.nf" }
process RUN_TRIPLE_PHENOTYPE_OVERLAP { script: "set -euo pipefail; ${params.nfbin} ${workflow.launchDir}/nf/triple_overlap.nf"}
process RUN_LD_CLUMPING { script: "set -euo pipefail; ${params.nfbin} ${workflow.launchDir}/nf/clump.nf"}
process RUN_COLOC { script: "set -euo pipefail; ${params.nfbin} ${workflow.launchDir}/nf/coloc.nf" }
process RUN_FUMA_PREP { script: "ste -euo pipefail; ${params.nfbin} ${workflow.launchDir}/nf/fuma_prep.nf"}

workflow {bvii
    RUN_QC()
    RUN_LDSC()
    RUN_LAVA()
    RUN_CONJFDR()
    RUN RUN_TRIPLE_PHENOTYPE_OVERLAP()
    RUN_LD_CLUMPING()
    RUN_COLOC()
    RUN_FUMA_PREP()
}
