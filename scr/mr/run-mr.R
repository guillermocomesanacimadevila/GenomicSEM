#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(remotes)
  library(progress)
  library(TwoSampleMR)
  library(phenoscanner)
  library(LDlinkR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript run_mr.R exposure.tsv outcome.tsv exposure_label outcome_label p_threshold ldlink_token confounder_file", call. = FALSE)
}

exp_file <- args[1]
out_file <- args[2]
exp_label <- args[3]
out_label <- args[4]
p_threshold <- as.numeric(args[5])
ld_token <- args[6]
conf_file <- args[7]

conf_terms <- scan(conf_file, what = "", quiet = TRUE)
if (length(conf_terms) == 0) stop("Confounder file is empty")
conf_pattern <- paste(conf_terms, collapse = "|")

outdir <- file.path(
  "../../outputs/MR", 
  paste0(
    exp_label, 
    "_", 
    out_label
    )
  )

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_exposure <- function(file, exposure_name) {
  dat <- read_exposure_data(
    filename = file,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FRQ",
    pval_col = "P",
    samplesize_col = "N",
    phenotype_col = "exposure",
    min_pval = 1e-200
  )
  dat$exposure <- exposure_name
  dat
}

read_outcome <- function(file, snps, outcome_name) {
  dat <- read_outcome_data(
    snps = snps,
    filename = file,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FRQ",
    pval_col = "P",
    samplesize_col = "N",
    phenotype_col = "outcome"
  )
  dat$outcome <- outcome_name
  dat
}

ld_clump_ldlink <- function(snps, token, pop = "EUR", r2 = 0.01) {
  keep <- c()
  remaining <- unique(snps)
  while (length(remaining) > 0) {
    lead <- remaining[1]
    keep <- c(keep, lead)
    ld <- tryCatch(LDproxy(snp = lead, pop = pop, r2d = "r2", token = token), error = function(e) NULL)
    if (is.null(ld)) {
      remaining <- setdiff(remaining, lead)
    } else {
      rm <- ld$RS_Number[ld$R2 >= r2]
      remaining <- setdiff(remaining, rm)
    }
  }
  keep
}

filter_confounders <- function(dat, pattern) {
  snps <- unique(dat$SNP)
  r <- tryCatch(phenoscanner(snp = snps), error = function(e) NULL)
  if (is.null(r) || is.null(r$results)) return(dat)
  hits <- r$results
  idx <- grep(pattern, hits$trait, ignore.case = TRUE)
  if (length(idx) == 0) return(dat)
  bad <- unique(hits$snp[idx])
  dat[!(dat$SNP %in% bad), ]
}

# confounders.txt e.g.
# smoking
# bmi
# obesity
# alcohol
# education
# socioeconomic
# blood pressure
# cholesterol
# triglyceride
# diabetes
# inflammation
# crp

run_one_direction <- function(exp_file, out_file, prefix, exp_name, out_name, p_threshold, outdir, token, pattern) {
  cat("Running", exp_name, "->", out_name, "\n")
  pb <- progress_bar$new(format = "[:bar] :percent :message", total = 6, clear = FALSE, width = 70)
  
  pb$tick(0, list(message = " loading exposure"))
  exp <- read_exposure(exp_file, exp_name)
  exp <- subset(exp, pval.exposure < p_threshold)
  if (nrow(exp) == 0) { pb$tick(6, list(message = " no instruments")); return(invisible()) }
  fwrite(as.data.table(exp), file.path(outdir, paste0(prefix, "_instruments_raw.csv")))
  
  exp$F <- (exp$beta.exposure^2) / (exp$se.exposure^2)
  exp$F_pval <- pchisq(exp$F, df = 1, lower.tail = FALSE)
  exp <- subset(exp, F > 10)
  if (nrow(exp) == 0) { pb$tick(6, list(message = " weak instruments")); return(invisible()) }
  pb$tick(1, list(message = " filtered instruments"))
  fwrite(as.data.table(exp[, c("SNP","beta.exposure","se.exposure","pval.exposure","F","F_pval")]),
         file.path(outdir, paste0(prefix, "_instruments_with_F.csv")))
  
  snp_keep <- ld_clump_ldlink(exp$SNP, token)
  exp <- exp[exp$SNP %in% snp_keep, ]
  if (nrow(exp) == 0) { pb$tick(6, list(message = " no LD-independent")); return(invisible()) }
  pb$tick(2, list(message = " LD clumped"))
  fwrite(as.data.table(exp), file.path(outdir, paste0(prefix, "_instruments_LDclumped.csv")))
  
  exp <- filter_confounders(exp, pattern)
  if (nrow(exp) == 0) { pb$tick(6, list(message = " all confounders")); return(invisible()) }
  pb$tick(3, list(message = " confounder-filtered"))
  fwrite(as.data.table(exp), file.path(outdir, paste0(prefix, "_instruments_no_confounders.csv")))
  
  out <- read_outcome(out_file, exp$SNP, out_name)
  harm <- harmonise_data(exp, out)
  if (nrow(harm) == 0) { pb$tick(6, list(message = " no harmonised")); return(invisible()) }
  pb$tick(4, list(message = " harmonised"))
  fwrite(as.data.table(harm), file.path(outdir, paste0(prefix, "_harmonised.csv")))
  
  mr_res <- mr(harm, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  het_res <- mr_heterogeneity(harm)
  pleio_res <- mr_pleiotropy_test(harm)
  loo_res <- mr_leaveoneout(harm)
  loo_plot <- mr_leaveoneout_plot(loo_res)
  
  fwrite(as.data.table(mr_res),    file.path(outdir, paste0(prefix, "_mr_results.csv")))
  fwrite(as.data.table(het_res),   file.path(outdir, paste0(prefix, "_heterogeneity.csv")))
  fwrite(as.data.table(pleio_res), file.path(outdir, paste0(prefix, "_pleiotropy.csv")))
  fwrite(as.data.table(loo_res),   file.path(outdir, paste0(prefix, "_leaveoneout.csv")))
  pdf(file.path(outdir, paste0(prefix, "_leaveoneout_plot.pdf"))); print(loo_plot[[1]]); dev.off()
  
  pb$tick(6, list(message = " done"))
}

run_one_direction(exp_file, out_file, "forward_expToOut", exp_label, out_label, p_threshold, outdir, ld_token, conf_pattern)
run_one_direction(out_file, exp_file, "reverse_outToExp", out_label, exp_label, p_threshold, outdir, ld_token, conf_pattern)
cat("Done\n")