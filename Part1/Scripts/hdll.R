#!/usr/bin/env Rscript

install.packages("remotes", repos = "https://cloud.r-project.org")
install.packages("RhpcBLASctl", repos = "https://cloud.r-project.org")
install.packages("progressr", repos = "https://cloud.r-project.org")

library(remotes)
packageVersion("remotes")

library(RhpcBLASctl)
packageVersion("RhpcBLASctl")

library(progressr)
handlers(global = TRUE)
handlers("txtprogressbar")

remotes::install_github("zhenin/HDL/HDL", dependencies = TRUE)
library(HDL)
library(dplyr)
library(data.table)
library(parallel)
packageVersion("HDL")

setwd("/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Checks")

# HDL-L per LD block 
LD.path <- "/Users/guillermocomesanacimadevila/LD_HDL_L/LD.path"
bim.path <- "/Users/guillermocomesanacimadevila/LD_HDL_L/bimfile"
stopifnot(file.exists(file.path(LD.path, "HDLL_LOC_snps.RData")))
load(file.path(LD.path, "HDLL_LOC_snps.RData")) 

sz_path <- "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Checks/SZ.hdll.tsv"
ad_path <- "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Checks/AD.hdll.tsv"

cov_min       <- 0.90       
write_snp_map <- TRUE        
out_dir       <- "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Checks/HDLL_per_block"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sz <- data.table::fread(sz_path)
ad <- data.table::fread(ad_path)
stopifnot(all(c("SNP","A1","A2","N","Z") %in% names(sz)))
stopifnot(all(c("SNP","A1","A2","N","Z") %in% names(ad)))

data.table::setDT(sz); data.table::setDT(ad)
sz[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]
ad[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]

data.table::setDT(NEWLOC)
NEWLOC[, `:=`(CHR = as.integer(CHR), piece = as.integer(piece))]
data.table::setkey(NEWLOC, CHR, piece)

block_bim_path <- function(chr, piece) {
  file.path(bim.path, sprintf("ukb_chr%d.%d_n336000.imputed_clean.bim", chr, piece))
}

block_info <- function(chr, piece, add_membership = FALSE) {
  bim_file <- block_bim_path(chr, piece)
  b <- data.table::fread(bim_file, header = FALSE, showProgress = FALSE)
  data.table::setnames(b, c("CHR","SNP","CM","BP","A1","A2"))
  n_ref <- nrow(b)
  present_sz <- b$SNP %chin% sz$SNP
  present_ad <- b$SNP %chin% ad$SNP
  cov_sz <- mean(present_sz)
  cov_ad <- mean(present_ad)
  bp_min <- if (n_ref) min(b$BP) else NA_integer_
  bp_max <- if (n_ref) max(b$BP) else NA_integer_
  
  cov_row <- data.table::data.table(
    CHR = chr, piece = piece, n_ref = n_ref,
    cov_scz = cov_sz, cov_ad = cov_ad,
    bp_min = bp_min, bp_max = bp_max
  )
  
  if (!add_membership) return(list(cov = cov_row, mem = NULL))
  
  mem <- data.table::data.table(
    CHR = chr, piece = piece, SNP = b$SNP,
    in_scz = as.integer(present_sz),
    in_ad  = as.integer(present_ad),
    BP = b$BP
  )
  list(cov = cov_row, mem = mem)
}

runner <- function(chr, piece) {
  HDL.L(
    gwas1 = sz, gwas2 = ad,
    Trait1name = "SCZ", Trait2name = "AD",
    LD.path = LD.path, bim.path = bim.path,
    chr = chr, piece = piece,
    eigen.cut = 0.99, N0 = 0, lim = exp(-18), output.file = ""
  )
}

cov_file      <- file.path(out_dir, "block_coverage.tsv")
snpmap_file   <- file.path(out_dir, "block_snp_membership.tsv.gz")
results_file  <- file.path(out_dir, "local_rg_per_block.tsv")
compact_file  <- file.path(out_dir, "local_rg_per_block_compact.csv")

with_progress({
  p <- progressr::progressor(steps = nrow(NEWLOC))
  cov_list <- vector("list", nrow(NEWLOC))
  mem_list <- if (write_snp_map) vector("list", nrow(NEWLOC)) else NULL
  
  for (i in seq_len(nrow(NEWLOC))) {
    chr_i   <- NEWLOC$CHR[i]
    piece_i <- NEWLOC$piece[i]
    p(sprintf("coverage chr %d piece %d", chr_i, piece_i))
    bi <- block_info(chr_i, piece_i, add_membership = write_snp_map)
    cov_list[[i]] <- bi$cov
    if (write_snp_map) mem_list[[i]] <- bi$mem
  }
  
  cov_tab <- data.table::rbindlist(cov_list)
  data.table::fwrite(cov_tab, cov_file, sep = "\t")
  
  if (write_snp_map) {
    mem_all <- data.table::rbindlist(mem_list, use.names = TRUE, fill = TRUE)
    data.table::fwrite(mem_all, snpmap_file, sep = "\t")
    rm(mem_all, mem_list); gc()
  }
})

message(sprintf("Coverage written: %s", cov_file))
if (write_snp_map) message(sprintf("SNP membership written: %s", snpmap_file))

cov_tab <- data.table::fread(cov_file)
idx_keep <- which(cov_tab$cov_scz >= cov_min & cov_tab$cov_ad >= cov_min)
if (!length(idx_keep)) {
  idx_keep <- seq_len(nrow(NEWLOC))
  message(sprintf("No blocks reached coverage ≥ %.0f%%. Proceeding with all (%d).",
                  100*cov_min, nrow(NEWLOC)))
} else {
  message(sprintf("Blocks kept by coverage ≥ %.0f%%: %d / %d",
                  100*cov_min, length(idx_keep), nrow(NEWLOC)))
}

with_progress({
  p2 <- progressr::progressor(steps = length(idx_keep))
  res_list <- lapply(idx_keep, function(i) {
    chr_i   <- NEWLOC$CHR[i]
    piece_i <- NEWLOC$piece[i]
    p2(sprintf("HDL-L chr %d piece %d", chr_i, piece_i))
    try(runner(chr_i, piece_i), silent = TRUE)
  })
  assign("res_list", res_list, inherits = TRUE)
})

ok <- vapply(res_list, function(x) !(inherits(x, "try-error") || is.null(x)), logical(1))
if (!any(ok)) stop("No blocks returned valid HDL-L results.")

res_local <- data.table::rbindlist(res_list[ok], fill = TRUE)

if (all(c("chr","piece") %in% names(res_local))) {
  data.table::setnames(res_local, c("chr","piece"), c("CHR","piece"))
}
res_local <- merge(res_local, cov_tab, by = c("CHR","piece"), all.x = TRUE)

if (all(c("Heritability_1","Heritability_2","Genetic_Covariance") %in% names(res_local))) {
  res_local[
    , local_rg := Genetic_Covariance / sqrt(pmax(Heritability_1, 0) * pmax(Heritability_2, 0))
  ]
} else {
  res_local[, local_rg := NA_real_]
}

data.table::setorder(res_local, CHR, piece)

data.table::fwrite(res_local, results_file, sep = "\t")

keep_cols <- intersect(c(
  "CHR","piece","bp_min","bp_max","n_ref","cov_scz","cov_ad",
  "Heritability_1","Heritability_2","Genetic_Covariance",
  "local_rg","eigen.cut","lim"
), names(res_local))
data.table::fwrite(res_local[, ..keep_cols], compact_file)

cat("\n=== DONE (HDL-L per LD block) ===\n")
cat(sprintf("Coverage table:      %s\n", cov_file))
if (write_snp_map) cat(sprintf("SNP membership:     %s\n", snpmap_file))
cat(sprintf("Per-block TSV:       %s\n", results_file))
cat(sprintf("Compact CSV:         %s\n", compact_file))
  


