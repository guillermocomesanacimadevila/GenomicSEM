#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(cfdr.pleio)
})

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
out_prefix <- args[2]

REF_DIR <- "/Users/c24102394/ref/conjFDR"
OUT_DIR <- "../../outputs/conjFDR"
LOCAL_REF_DIR <- paste0(out_prefix, "_localref")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- fread(infile)

trait_scz <- dat[, .(SNP, BETA = BETA_SCZ, PVAL = P_SCZ)]
trait_ad  <- dat[, .(SNP, BETA = BETA_AD, PVAL = P_AD)]

obj <- cfdr_pleio$new()

obj$init_data(
  trait1 = trait_scz,
  trait2 = trait_ad,
  trait_names = c("SCZ", "AD"),
  refdat = refdata_location(REF_DIR),
  local_refdat_path = LOCAL_REF_DIR,
  verbose = FALSE
)

obj$initialize_pruning_index(n_iter = 50, seed = 154226, verbose = FALSE)
obj$calculate_cond_fdr(fdr_trait = 1, verbose = FALSE)
obj$calculate_cond_fdr(fdr_trait = 2, verbose = FALSE)

res <- obj$get_trait_results()

fwrite(res, file.path(OUT_DIR, paste0(out_prefix, "_cfdr_results.tsv")), sep = "\t")
fwrite(res[conj_fdr < 0.05], file.path(OUT_DIR, paste0(out_prefix, "_shared_hits.tsv")), sep = "\t")
