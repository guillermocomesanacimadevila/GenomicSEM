#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(progress)
  library(data.table)
  library(remotes)
  library(phenoscanner)
  library(TwoSampleMR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript run-mr.R exposure.tsv outcome.tsv exposure_label outcome_label p_threshold ldlink_token confounder_file",
       call. = False
  )
}

exp_file <- args[1]
out_file <- args[2]
exp_label <- args[3]
out_label <- args[4]
p_threshold <- as.numeric(args[5])
ld_token <- args[6]
conf_file <- args[7]

outdir <- file.path("../../outputs/MR/",
                    paste0(exp_label, "_", out_label
                           )
                    )
dir.create(
  outdir, 
  showWarnings = FALSE,
  recursive = TRUE
  )

read_exposure <- function(file, exposure_name) {
  d <- suppressMessages(
    read_exposure_data(
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
  )
  d$exposure <- exposure_name
  d
}

read_outcome <- function(file, outcome_name) {
  d <- suppressMessages(
    read_outcome_data(
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
      phenotype_col = "outcome",
      min_pval = 1e-200
    )
  )
  d$outcome <- outcome_name
  d
}

ld_clump_ldlink <- function(dat, clump_kb = 1000, clump_r2 = 0.01) {
  total <- nrow(dat)
  cat(sprintf("[3/6] LD clumping: %d SNPs before clumping\n", total))
  res <- clump_data(
    dat,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p1 = 1,
    clump_p2 = 1
  )
  kept <- nrow(res)
  cat(sprintf("[3/6] LD clumping: %d SNPs after clumping\n", kept))
  res
}


filter_confounders <- function(dat, pattern) {
  snps <- unique(dat$SNP)
  if (length(snps) == 0) return(dat)
  
  chunk_size <- 10L
  n_chunks <- ceiling(length(snps) / chunk_size)
  res_list <- vector("list", n_chunks)
  
  for (i in seq_len(n_chunks)) {
    idx <- ((i - 1L) * chunk_size + 1L):min(i * chunk_size, length(snps))
    this_snps <- snps[idx]
    
    cat(sprintf("\r[4/6] Phenoscanner confounder filtering: batch %d/%d", i, n_chunks))
    flush.console()
    
    r <- tryCatch(
      phenoscanner(snp = this_snps),
      error = function(e) NULL
    )
    
    if (!is.null(r) && !is.null(r$results)) {
      res_list[[i]] <- r$results
    }
  }
  
  cat("\n")  
  
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (length(res_list) == 0) return(dat)
  
  hits <- do.call(rbind, res_list)
  idx  <- grep(pattern, hits$trait, ignore.case = TRUE)
  if (length(idx) == 0) return(dat)
  
  bad <- unique(hits$snp[idx])
  dat[!(dat$SNP %in% bad), ]
}


run_one_direction <- function(exp_file, out_file, prefix, exp_name, out_name, p_threshold, outdir, token, pattern) {
  cat("=== Running", exp_name, "->", out_name, "===\n")
  
  cat("[1/6] Loading exposure and p-value filtering\n")
  exp <- read_exposure(exp_file, exp_name)
  exp <- subset(exp, pval.exposure < p_threshold)
  if (nrow(exp) == 0) {
    cat("No instruments at p <", p_threshold, "for", exp_name, "\n")
    return(invisible())
  }
  fwrite(as.data.table(exp),
         file.path(outdir, paste0(prefix, "_instruments_raw.csv")))
  
  cat("[2/6] F-stat calculation and weak instrument filter\n")
  exp$F <- (exp$beta.exposure^2) / (exp$se.exposure^2)
  exp$F_pval <- pchisq(exp$F, df = 1, lower.tail = FALSE)
  exp <- subset(exp, F > 10)
  if (nrow(exp) == 0) {
    cat("All instruments for", exp_name, "are weak (F <= 10)\n")
    return(invisible())
  }
  fwrite(as.data.table(exp[, c("SNP","beta.exposure","se.exposure","pval.exposure","F","F_pval")]),
         file.path(outdir, paste0(prefix, "_instruments_with_F.csv")))
  
  cat("[3/6] LD clumping\n")
  exp <- ld_clump_ldlink(exp)
  if (nrow(exp) == 0) {
    cat("No LD-independent instruments remain for", exp_name, "\n")
    return(invisible())
  }
  fwrite(as.data.table(exp),
         file.path(outdir, paste0(prefix, "_instruments_LDclumped.csv")))
  
  cat("[4/6] Phenoscanner confounder filtering\n")
  exp <- filter_confounders(exp, pattern)
  if (nrow(exp) == 0) {
    cat("All LD-independent instruments flagged as confounder-linked for", exp_name, "\n")
    return(invisible())
  }
  fwrite(as.data.table(exp),
         file.path(outdir, paste0(prefix, "_instruments_no_confounders.csv")))
  
  cat("[5/6] Harmonising exposure and outcome\n")
  out <- read_outcome(out_file, out_name)
  out <- out[out$SNP %in% exp$SNP, ]
  harm <- harmonise_data(exp, out)
  if (nrow(harm) == 0) {
    cat("No harmonised SNPs for", exp_name, "->", out_name, "\n")
    return(invisible())
  }
  fwrite(as.data.table(harm),
         file.path(outdir, paste0(prefix, "_harmonised.csv")))
  
  cat("[6/6] Running MR models and sensitivity analyses\n")
  mr_res <- mr(harm, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  mr_res$lo_ci <- mr_res$b - 1.96 * mr_res$se
  mr_res$up_ci <- mr_res$b + 1.96 * mr_res$se
  mr_res$OR    <- exp(mr_res$b)
  mr_res$OR_lo <- exp(mr_res$lo_ci)
  mr_res$OR_up <- exp(mr_res$up_ci)
  
  het_res <- mr_heterogeneity(harm)
  pleio_res <- mr_pleiotropy_test(harm)
  loo_res <- mr_leaveoneout(harm)
  loo_plot <- mr_leaveoneout_plot(loo_res)
  
  fwrite(as.data.table(mr_res),
         file.path(outdir, paste0(prefix, "_mr_results.csv")))
  fwrite(as.data.table(het_res),
         file.path(outdir, paste0(prefix, "_heterogeneity.csv")))
  fwrite(as.data.table(pleio_res),
         file.path(outdir, paste0(prefix, "_pleiotropy.csv")))
  fwrite(as.data.table(loo_res),
         file.path(outdir, paste0(prefix, "_leaveoneout.csv")))
  
  pdf(file.path(outdir, paste0(prefix, "_leaveoneout_plot.pdf")))
  print(loo_plot[[1]])
  dev.off()
  
  cat("=== Finished", exp_name, "->", out_name, "===\n")
}

run_one_direction(exp_file, out_file, "forward_expToOut", exp_label, out_label, p_threshold, outdir, ld_token, conf_pattern)
run_one_direction(out_file, exp_file, "reverse_outToExp", out_label, exp_label, p_threshold, outdir, ld_token, conf_pattern)

cat("Done\n")