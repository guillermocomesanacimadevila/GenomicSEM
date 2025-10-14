#!/usr/bin/env RScript
# Check rg between AD and SZ

install.packages(c("remotes", "gitcreds"))

library(gitcreds)
try(gitcreds_delete(), silent = TRUE)  

Sys.unsetenv("GITHUB_PAT")
Sys.unsetenv("GITHUB_TOKEN")
options(github_pat = NULL)

# Install GenomicSEM without a GitHub token
library(remotes)
remotes::install_github("GenomicSEM/GenomicSEM", auth_token = NULL, upgrade = "never")

library(GenomicSEM)
packageVersion("GenomicSEM")

install.packages("RhpcBLASctl")
library(RhpcBLASctl)
packageVersion("RhpcBLASctl")

install.packages("progressr")
library(progressr)
handlers(global = TRUE)
handlers("txtprogressbar")

# HDL-L package - check for genetic correlation with this one vs LDSC
remotes::install_github("zhenin/HDL/HDL", dependencies = TRUE)

library(HDL)
library(dplyr)
library(data.table)
library(parallel)
packageVersion("HDL")

# You already did this earlier, but do it again to be explicit:
setwd("/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data")

# 1) Run munge (as you have)
munged <- munge(
  files       = c("SZ/SZ/SZ.ldsc.sumstats", "AD/AD/AD.ldsc.sumstats"),
  hm3         = "/Users/guillermocomesanacimadevila/ldsc_ref/w_hm3.snplist",
  trait.names = c("SCZ","AD"),
  N           = c(171880, 57693),
  info.filter = 0.9, 
  maf.filter  = 0.01  
)

# 2) Construct absolute paths to the munged files (they’re saved in getwd())
out_files <- file.path(getwd(), c("SCZ.sumstats.gz", "AD.sumstats.gz"))

# 3) Sanity check they exist BEFORE calling ldsc()
print(out_files)
stopifnot(all(file.exists(out_files)))  # will error if not found

# 4) Point to LD score resources
ld_path  <- "/Users/guillermocomesanacimadevila/ldsc_ref/eur_w_ld_chr"
wld_path <- ld_path

# 5) Run LDSC
# sample.prev = cases / (cases + controls)
ldsc_out <- ldsc(
  traits          = out_files,
  sample.prev     = c(0.5, 0.5),
  population.prev = c(0.01, 0.05),
  ld              = ld_path,
  wld             = wld_path,
  trait.names     = c("SCZ","AD")
)

ldsc_out

# Save these 
write.csv(ldsc_out$V, "ldsc_V.csv", row.names = FALSE)
write.csv(ldsc_out$S, "ldsc_S.csv", row.names = FALSE)
write.csv(ldsc_out$I, "ldsc_I.csv", row.names = FALSE)
write.csv(ldsc_out$N, "ldsc_N.csv", row.names = FALSE)
write.csv(data.frame(m = ldsc_out$m), "ldsc_m.csv", row.names = FALSE)

ldsc_rg_summary <- function(ldsc_out, conf.level = 0.95) {
  S <- ldsc_out$S; V <- ldsc_out$V
  rg <- S[1,2] / sqrt(S[1,1] * S[2,2])
  g  <- c(-0.5*rg/S[1,1], 1/sqrt(S[1,1]*S[2,2]), -0.5*rg/S[2,2])
  se <- sqrt(as.numeric(t(g) %*% V %*% g))
  z  <- rg / se
  p  <- 2 * pnorm(-abs(z))
  zcrit <- qnorm(1 - (1 - conf.level)/2)
  ci <- c(max(-1, rg - zcrit*se), min(1, rg + zcrit*se))
  data.frame(rg = rg, SE = se, z = z, p = p, CI_low = ci[1], CI_high = ci[2])
}

ldsc_rg_summary(ldsc_out)

# Run 1 -> rg = 0.08694406, SE = 0.04080815, z = 2.131, p = 0.03313
# Run 1 -> Neff (Not Ntotal) - change sample.prev accordingly 

# Run 2 -> rg = 0.08694406, SE = 0.04080815, z = 2.131, p = 0.03313
# Run 1 -> Ntotal - changed sample.prev  

# Run 3 -> rg = 0.08355765, SE = 0.04130047, z = 2.023, p = 0.04306
# Run 3 -> Proper INDEL removal from AD GWAS, Apply INFO > 0.75 on SZ & MAF < 0.01 ON BOTH

# Run 4 -> rg = 0.08258438, SE = 0.04301183, z = 1.920, p = 0.05485
# Run 4 -> Change INFO cut-off from 0.75 to 0.90 

# Run 5 -> rg = 0.08258438, SE = 0.04301183, z = 1.920, p = 0.05485, CI (95%) = -0.00171726 +/- 0.166886
# Run 5 -> Back to Neff


# ===================== # 
# HDL-L -> recompute rg

# 1. Set LD_file and bim file 
LD.path <- "/Users/guillermocomesanacimadevila/LD_HDL_L/LD.path"
bim.path <- "/Users/guillermocomesanacimadevila/LD_HDL_L/bimfile"
stopifnot(file.exists(file.path(LD.path, "HDLL_LOC_snps.RData")))
load(file.path(LD.path, "HDLL_LOC_snps.RData")) 

# Re-harmonised GWAS summary stats for HDL-L in .py 
sz_gwas <- read.csv("SZ/SZ.hdll.tsv", sep="\t")
ad_gwas <- read.csv("AD/AD.hdll.tsv", sep="\t")

print(dim(sz_gwas))
print(dim(ad_gwas))

# Sanity cehcks 
stopifnot(all(c("SNP","A1","A2","N","Z") %in% names(sz_gwas)))
stopifnot(all(c("SNP","A1","A2","N","Z") %in% names(ad_gwas)))
sz_gwas$N <- as.numeric(sz_gwas$N);  sz_gwas$Z <- as.numeric(sz_gwas$Z)
ad_gwas$N <- as.numeric(ad_gwas$N);  ad_gwas$Z <- as.numeric(ad_gwas$Z)


# Test run 
chr_test   <- as.numeric(NEWLOC$CHR[1])
piece_test <- as.numeric(NEWLOC$piece[1])

res_test <- HDL.L(
  gwas1=sz_gwas, gwas2=ad_gwas,
  Trait1name="SCZ", Trait2name="AD",
  LD.path=LD.path, bim.path=bim.path,
  chr=chr_test, piece=piece_test,
  eigen.cut=0.99, N0=0, lim=exp(-18), output.file=""
)

res_df <- as.data.frame(res_test)
View(res_df)

## ========== WHOLE GENOME NOW (no manual core settings) ========== #

# ensure DT + numeric
setDT(sz_gwas); setDT(ad_gwas)
sz_gwas[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]
ad_gwas[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]

# ensure NEWLOC integer
setDT(NEWLOC)
NEWLOC[, `:=`(CHR = as.integer(CHR), piece = as.integer(piece))]

# coverage helper
locus_coverage <- function(chr, piece) {
  bim_file <- file.path(bim.path, sprintf("ukb_chr%d.%d_n336000.imputed_clean.bim", chr, piece))
  b <- data.table::fread(bim_file, header = FALSE, showProgress = FALSE)
  data.table::setnames(b, c("CHR","SNP","CM","BP","A1","A2"))
  n_ref   <- nrow(b)
  cov_scz <- mean(b$SNP %chin% sz_gwas$SNP)
  cov_ad  <- mean(b$SNP %chin% ad_gwas$SNP)
  data.table(CHR = chr, piece = piece, n_ref = n_ref, cov_scz = cov_scz, cov_ad = cov_ad)
}

# 1) pre-compute coverage (sequential, with progress)
with_progress({
  p <- progressor(steps = nrow(NEWLOC))
  cov_list <- lapply(seq_len(nrow(NEWLOC)), function(i) {
    p(sprintf("coverage chr %d piece %d", NEWLOC$CHR[i], NEWLOC$piece[i]))
    locus_coverage(NEWLOC$CHR[i], NEWLOC$piece[i])
  })
  assign("cov_list", cov_list, inherits = TRUE)
})
cov_tab <- rbindlist(cov_list)

# 2) pick loci by coverage
cov_min  <- 0.90
idx_all  <- seq_len(nrow(NEWLOC))
idx_keep <- idx_all[(cov_tab$cov_scz >= cov_min) & (cov_tab$cov_ad >= cov_min)]
if (length(idx_keep) == 0L) idx_keep <- idx_all
message(sprintf("Loci kept by coverage ≥ %.0f%%: %d / %d",
                100*cov_min, length(idx_keep), nrow(NEWLOC)))

# 3) HDL-L runner
runner <- function(chr, piece) {
  HDL.L(
    gwas1 = sz_gwas, gwas2 = ad_gwas,
    Trait1name = "SCZ", Trait2name = "AD",
    LD.path = LD.path, bim.path = bim.path,
    chr = chr, piece = piece,
    eigen.cut = 0.99, N0 = 0, lim = exp(-18), output.file = ""
  )
}

# 4) run HDL-L over kept loci (sequential, with progress)
with_progress({
  p2 <- progressor(steps = length(idx_keep))
  res_list <- lapply(idx_keep, function(i) {
    chr   <- NEWLOC$CHR[i]
    piece <- NEWLOC$piece[i]
    p2(sprintf("HDL-L chr %d piece %d", chr, piece))
    try(runner(chr, piece), silent = TRUE)
  })
  assign("res_list", res_list, inherits = TRUE)
})

# 5) collect & merge with coverage
ok <- vapply(res_list, function(x) !(inherits(x, "try-error") || is.null(x)), logical(1))
if (!any(ok)) {
  dir.create("HDLL_out", showWarnings = FALSE)
  fwrite(cov_tab, "HDLL_out/SCZ_AD.HDLL.coverage_only.tsv", sep = "\t")
  stop("No loci returned valid HDL-L results. Wrote coverage table for debugging.")
}
res_local <- rbindlist(res_list[ok], fill = TRUE)

if (all(c("chr","piece") %in% names(res_local))) {
  setnames(res_local, c("chr","piece"), c("CHR","piece"))
}
res_join <- merge(res_local, cov_tab, by = c("CHR","piece"), all.x = TRUE)
if (all(c("CHR","piece") %in% names(res_join))) {
  setnames(res_join, c("CHR","piece"), c("chr","piece"))
}

# 6) save local results
dir.create("HDLL_out", showWarnings = FALSE)
fwrite(res_join, "HDLL_out/SCZ_AD.HDLL.local_rg.with_coverage.tsv", sep = "\t")
saveRDS(res_join, "HDLL_out/SCZ_AD.HDLL.local_rg.with_coverage.rds")
print(dim(res_join))

# 7) global rg on valid, covered loci
valid <- res_join[
  is.finite(Heritability_1) & is.finite(Heritability_2) & is.finite(Genetic_Covariance) &
    Heritability_1 > 0 & Heritability_2 > 0 &
    cov_scz >= cov_min & cov_ad >= cov_min
]

if (nrow(valid) > 1L) {
  h1   <- valid$Heritability_1
  h2   <- valid$Heritability_2
  gcov <- valid$Genetic_Covariance
  rg_hat <- sum(gcov) / sqrt(sum(h1) * sum(h2))
  
  theta_i <- sapply(seq_len(nrow(valid)), function(i) {
    h1_i <- sum(h1) - h1[i]
    h2_i <- sum(h2) - h2[i]
    gcov_i <- sum(gcov) - gcov[i]
    if (h1_i <= 0 || h2_i <= 0) return(NA_real_)
    gcov_i / sqrt(h1_i * h2_i)
  })
  theta_i <- theta_i[is.finite(theta_i)]
  n <- length(theta_i)
  theta_dot <- mean(theta_i)
  se_jk <- sqrt(((n - 1) / n) * sum((theta_i - theta_dot)^2))
  z <- rg_hat / se_jk
  p <- 2 * pnorm(-abs(z))
  ci_low  <- max(-1, rg_hat - 1.96 * se_jk)
  ci_high <- min(1, rg_hat + 1.96 * se_jk)
  
  global_tab <- data.table(
    rg = rg_hat, SE = se_jk, z = z, p = p,
    CI_low = ci_low, CI_high = ci_high,
    n_loci_used = nrow(valid),
    cov_min_used = cov_min
  )
  fwrite(global_tab, "HDLL_out/SCZ_AD.HDLL.global_rg.tsv", sep = "\t")
  saveRDS(global_tab, "HDLL_out/SCZ_AD.HDLL.global_rg.rds")
  print(global_tab)
} else {
  message("Not enough valid loci (positive h2 & finite gcov under coverage rule) to compute a global rg.")
}

# ========== no mhc ========== #
## Re-run LDSC without MHC complex - remove from both GWAS 
# Drop that in py (already created: SZ/SZ.noMHC.ldsc.sumstats, AD/AD.noMHC.ldsc.sumstats)

sz_noMHC <- "SZ/SZ.ldsc.clean.noMHC.sumstats"
ad_noMHC <- "AD/AD.ldsc.clean.noMHC.sumstats"
stopifnot(file.exists(sz_noMHC), file.exists(ad_noMHC))

count_rows <- function(path) data.table::fread(path, showProgress = FALSE)[, .N]
n_in_sz <- count_rows(sz_noMHC)
n_in_ad <- count_rows(ad_noMHC)

munged_noMHC <- munge(
  files       = c(sz_noMHC, ad_noMHC),
  hm3         = "/Users/guillermocomesanacimadevila/ldsc_ref/w_hm3.snplist",
  trait.names = c("SCZ_noMHC","AD_noMHC"),
  N           = c(171880, 57693),
  info.filter = 0.90,
  maf.filter  = 0.01
)

out_files_noMHC <- file.path(getwd(), c("SCZ_noMHC.sumstats.gz", "AD_noMHC.sumstats.gz"))
stopifnot(all(file.exists(out_files_noMHC)))
n_mg_sz <- count_rows(out_files_noMHC[1])
n_mg_ad <- count_rows(out_files_noMHC[2])

ldsc_noMHC <- ldsc(
  traits          = out_files_noMHC,
  sample.prev     = c(0.5, 0.5),
  population.prev = c(0.01, 0.05),
  ld              = ld_path,
  wld             = wld_path,
  trait.names     = c("SCZ_noMHC","AD_noMHC")
)

print(ldsc_noMHC)
print(ldsc_rg_summary(ldsc_noMHC))

cat(sprintf("\nSNPs used by LDSC (m): %,d\n\n", ldsc_noMHC$m))

write.csv(ldsc_noMHC$V, "ldsc_noMHC_V.csv", row.names = FALSE)
write.csv(ldsc_noMHC$S, "ldsc_noMHC_S.csv", row.names = FALSE)
write.csv(ldsc_noMHC$I, "ldsc_noMHC_I.csv", row.names = FALSE)
write.csv(ldsc_noMHC$N, "ldsc_noMHC_N.csv", row.names = FALSE)
write.csv(data.frame(m = ldsc_noMHC$m), "ldsc_noMHC_m.csv", row.names = FALSE)
write.csv(ldsc_rg_summary(ldsc_noMHC), "ldsc_noMHC_rg_summary.csv", row.names = FALSE)

counts_tbl <- data.frame(
  Stage = c("Input_noMHC","Input_noMHC","PostMunge","PostMunge","LDSC_m"),
  Trait = c("SCZ","AD","SCZ","AD","Both"),
  SNPs  = c(n_in_sz, n_in_ad, n_mg_sz, n_mg_ad, ldsc_noMHC$m)
)
print(counts_tbl)
write.csv(counts_tbl, "noMHC_snp_counts_summary.csv", row.names = FALSE)

# ===== 
# LDSC without chromosome 6 
sz_noChr6 <- "SZ/SZ/SZ.ldsc.clean.noChr6.sumstats"
ad_noChr6 <- "AD/AD/AD.ldsc.clean.noChr6.sumstats"
stopifnot(file.exists(sz_noChr6), file.exists(ad_noChr6))

n_in_sz <- count_rows(sz_noChr6)
n_in_ad <- count_rows(ad_noChr6)

munged_noChr6 <- munge(
  files       = c(sz_noChr6, ad_noChr6),
  hm3         = "/Users/guillermocomesanacimadevila/ldsc_ref/w_hm3.snplist",
  trait.names = c("SCZ_noChr6","AD_noChr6"),
  N           = c(171880, 57693),
  info.filter = 0.90,
  maf.filter  = 0.01
)

out_files_noChr6 <- file.path(getwd(), c("SCZ_noChr6.sumstats.gz", "AD_noChr6.sumstats.gz"))
stopifnot(all(file.exists(out_files_noChr6)))
n_mg_sz <- count_rows(out_files_noChr6[1])
n_mg_ad <- count_rows(out_files_noChr6[2])

cat("\n[noChr6 SNP counts]\n")
cat(sprintf("SCZ: input = %s  -> post-munge = %s\n",
            format(n_in_sz, big.mark=","), format(n_mg_sz, big.mark=",")))
cat(sprintf("AD : input = %s  -> post-munge = %s\n\n",
            format(n_in_ad, big.mark=","), format(n_mg_ad, big.mark=",")))

ldsc_noChr6 <- ldsc(
  traits          = out_files_noChr6,
  sample.prev     = c(0.5, 0.5),
  population.prev = c(0.01, 0.05),
  ld              = ld_path,
  wld             = wld_path,
  trait.names     = c("SCZ_noChr6","AD_noChr6")
)

print(ldsc_noChr6)
print(ldsc_rg_summary(ldsc_noChr6))

cat(sprintf("\nSNPs used by LDSC (m): %s\n\n", format(ldsc_noChr6$m, big.mark=",")))

write.csv(ldsc_noChr6$V, "ldsc_noChr6_V.csv", row.names = FALSE)
write.csv(ldsc_noChr6$S, "ldsc_noChr6_S.csv", row.names = FALSE)
write.csv(ldsc_noChr6$I, "ldsc_noChr6_I.csv", row.names = FALSE)
write.csv(ldsc_noChr6$N, "ldsc_noChr6_N.csv", row.names = FALSE)
write.csv(data.frame(m = ldsc_noChr6$m), "ldsc_noChr6_m.csv", row.names = FALSE)
write.csv(ldsc_rg_summary(ldsc_noChr6), "ldsc_noChr6_rg_summary.csv", row.names = FALSE)

counts_tbl <- data.frame(
  Stage = c("Input_noChr6","Input_noChr6","PostMunge","PostMunge","LDSC_m"),
  Trait = c("SCZ","AD","SCZ","AD","Both"),
  SNPs  = c(n_in_sz, n_in_ad, n_mg_sz, n_mg_ad, ldsc_noChr6$m)
)
print(counts_tbl)
write.csv(counts_tbl, "noChr6_snp_counts_summary.csv", row.names = FALSE)


# Run HDL-L with no-MHC
# Prepare for HDL-L in .py
sz_noMHC_hdl <- read.csv("SZ/SZ.hdl.sumstats", sep="\t")
ad_noMHC_hdl <- read.csv("AD/AD.hdl.sumstats", sep="\t")
print(dim(sz_noMHC_hdl))
print(dim(ad_noMHC_hdl))

stopifnot(all(c("SNP","A1","A2","N","Z") %in% names(sz_noMHC_hdl)))
stopifnot(all(c("SNP","A1","A2","N","Z") %in% names(ad_noMHC_hdl)))
sz_noMHC_hdl$N <- as.numeric(sz_noMHC_hdl$N);  sz_noMHC_hdl$Z <- as.numeric(sz_noMHC_hdl$Z)
ad_noMHC_hdl$N <- as.numeric(ad_noMHC_hdl$N);  ad_noMHC_hdl$Z <- as.numeric(ad_noMHC_hdl$Z)

chr_test_noMHC <- as.numeric(NEWLOC$CHR[1])
piece_test_noMHC <- as.numeric(NEWLOC$piece[1])

res_test_noMHC <- HDL.L(
  gwas1=sz_noMHC_hdl, gwas2=ad_noMHC_hdl,
  Trait1name="SCZ", Trait2name="AD",
  LD.path=LD.path, bim.path=bim.path,
  chr=chr_test_noMHC, piece=piece_test_noMHC,
  eigen.cut=0.99, N0=0, lim=exp(-18), output.file=""
)

res_df_noMHC <- as.data.frame(res_test_noMHC)
View(res_df_noMHC)

setDT(sz_noMHC_hdl); setDT(ad_noMHC_hdl)
sz_noMHC_hdl[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]
ad_noMHC_hdl[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]

setDT(NEWLOC)
NEWLOC[, `:=`(CHR = as.integer(CHR), piece = as.integer(piece))]

## ===== HDL-L (no-MHC): coverage -> per-locus -> aggregate -> global rg =====
locus_coverage_noMHC <- function(chr, piece) {
  bim_file <- file.path(bim.path, sprintf("ukb_chr%d.%d_n336000.imputed_clean.bim", chr, piece))
  b <- data.table::fread(bim_file, header = FALSE, showProgress = FALSE)
  data.table::setnames(b, c("CHR","SNP","CM","BP","A1","A2"))
  n_ref   <- nrow(b)
  cov_scz <- mean(b$SNP %chin% sz_noMHC_hdl$SNP)
  cov_ad  <- mean(b$SNP %chin% ad_noMHC_hdl$SNP)
  data.table(CHR = chr, piece = piece, n_ref = n_ref, cov_scz = cov_scz, cov_ad = cov_ad)
}

with_progress({
  p <- progressor(steps = nrow(NEWLOC))
  cov_list_noMHC <- lapply(seq_len(nrow(NEWLOC)), function(i) {
    p(sprintf("coverage noMHC chr %d piece %d", NEWLOC$CHR[i], NEWLOC$piece[i]))
    locus_coverage_noMHC(NEWLOC$CHR[i], NEWLOC$piece[i])
  })
})
cov_tab_noMHC <- rbindlist(cov_list_noMHC)

cov_min_noMHC <- 0.90
idx_all  <- seq_len(nrow(NEWLOC))
idx_keep <- idx_all[(cov_tab_noMHC$cov_scz >= cov_min_noMHC) & (cov_tab_noMHC$cov_ad >= cov_min_noMHC)]
if (length(idx_keep) == 0L) idx_keep <- idx_all
message(sprintf("[noMHC] Loci kept by coverage ≥ %.0f%%: %d / %d",
                100*cov_min_noMHC, length(idx_keep), nrow(NEWLOC)))

runner_noMHC <- function(chr, piece) {
  HDL.L(
    gwas1 = sz_noMHC_hdl, gwas2 = ad_noMHC_hdl,
    Trait1name = "SCZ", Trait2name = "AD",
    LD.path = LD.path, bim.path = bim.path,
    chr = chr, piece = piece,
    eigen.cut = 0.99, N0 = 0, lim = exp(-18), output.file = ""
  )
}

with_progress({
  p2 <- progressor(steps = length(idx_keep))
  res_list_noMHC <- lapply(idx_keep, function(i) {
    chr   <- NEWLOC$CHR[i]
    piece <- NEWLOC$piece[i]
    p2(sprintf("HDL-L noMHC chr %d piece %d", chr, piece))
    try(runner_noMHC(chr, piece), silent = TRUE)
  })
})

ok_noMHC <- vapply(res_list_noMHC, function(x) !(inherits(x, "try-error") || is.null(x)), logical(1))
dir.create("HDLL_out_noMHC", showWarnings = FALSE)
if (!any(ok_noMHC)) {
  fwrite(cov_tab_noMHC, "HDLL_out_noMHC/SCZ_AD.noMHC.HDLL.coverage_only.tsv", sep = "\t")
  stop("[noMHC] No loci returned valid HDL-L results. Wrote coverage table for debugging.")
}
res_local_noMHC <- rbindlist(res_list_noMHC[ok_noMHC], fill = TRUE)

if (all(c("chr","piece") %in% names(res_local_noMHC))) {
  setnames(res_local_noMHC, c("chr","piece"), c("CHR","piece"))
}
res_join_noMHC <- merge(res_local_noMHC, cov_tab_noMHC, by = c("CHR","piece"), all.x = TRUE)
if (all(c("CHR","piece") %in% names(res_join_noMHC))) {
  setnames(res_join_noMHC, c("CHR","piece"), c("chr","piece"))
}

fwrite(res_join_noMHC, "HDLL_out_noMHC/SCZ_AD.noMHC.HDLL.local_rg.with_coverage.tsv", sep = "\t")
saveRDS(res_join_noMHC, "HDLL_out_noMHC/SCZ_AD.noMHC.HDLL.local_rg.with_coverage.rds")
print(dim(res_join_noMHC))

valid_noMHC <- res_join_noMHC[
  is.finite(Heritability_1) & is.finite(Heritability_2) & is.finite(Genetic_Covariance) &
    Heritability_1 > 0 & Heritability_2 > 0 &
    cov_scz >= cov_min_noMHC & cov_ad >= cov_min_noMHC
]

if (nrow(valid_noMHC) > 1L) {
  h1   <- valid_noMHC$Heritability_1
  h2   <- valid_noMHC$Heritability_2
  gcov <- valid_noMHC$Genetic_Covariance
  
  rg_hat <- sum(gcov) / sqrt(sum(h1) * sum(h2))
  
  theta_i <- sapply(seq_len(nrow(valid_noMHC)), function(i) {
    h1_i   <- sum(h1) - h1[i]
    h2_i   <- sum(h2) - h2[i]
    gcov_i <- sum(gcov) - gcov[i]
    if (h1_i <= 0 || h2_i <= 0) return(NA_real_)
    gcov_i / sqrt(h1_i * h2_i)
  })
  theta_i <- theta_i[is.finite(theta_i)]
  n <- length(theta_i)
  theta_dot <- mean(theta_i)
  se_jk <- sqrt(((n - 1) / n) * sum((theta_i - theta_dot)^2))
  z <- rg_hat / se_jk
  p <- 2 * pnorm(-abs(z))
  ci_low  <- max(-1, rg_hat - 1.96 * se_jk)
  ci_high <- min(1, rg_hat + 1.96 * se_jk)
  
  global_tab_noMHC <- data.table(
    rg = rg_hat, SE = se_jk, z = z, p = p,
    CI_low = ci_low, CI_high = ci_high,
    n_loci_used = nrow(valid_noMHC),
    cov_min_used = cov_min_noMHC
  )
  fwrite(global_tab_noMHC, "HDLL_out_noMHC/SCZ_AD.noMHC.HDLL.global_rg.tsv", sep = "\t")
  saveRDS(global_tab_noMHC, "HDLL_out_noMHC/SCZ_AD.noMHC.HDLL.global_rg.rds")
  print(global_tab_noMHC)
} else {
  message("[noMHC] Not enough valid loci (positive h2 & finite gcov under coverage rule) to compute a global rg.")
}




# ===== NOW WITHOUT CHROMOSOME 6 ===== #
sz_nochr6_hdl <- read.csv("SZ/SZ.nochr6hdl.sumstats", sep="\t")
ad_nochr6_hdl <- read.csv("AD/AD.nochr6hdl.sumstats", sep="\t")
print(dim(sz_nochr6_hdl))
print(dim(ad_nochr6_hdl))

stopifnot(all(c("SNP","A1","A2","N","Z") %in% names(sz_nochr6_hdl)))
stopifnot(all(c("SNP","A1","A2","N","Z") %in% names(ad_nochr6_hdl)))
setDT(sz_nochr6_hdl); setDT(ad_nochr6_hdl)
sz_nochr6_hdl[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]
ad_nochr6_hdl[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]

setDT(NEWLOC)
NEWLOC[, `:=`(CHR = as.integer(CHR), piece = as.integer(piece))]
NEWLOC_no6 <- NEWLOC[CHR != 6L]

locus_coverage_noChr6 <- function(chr, piece) {
  bim_file <- file.path(bim.path, sprintf("ukb_chr%d.%d_n336000.imputed_clean.bim", chr, piece))
  b <- data.table::fread(bim_file, header = FALSE, showProgress = FALSE)
  data.table::setnames(b, c("CHR","SNP","CM","BP","A1","A2"))
  n_ref   <- nrow(b)
  cov_scz <- mean(b$SNP %chin% sz_nochr6_hdl$SNP)
  cov_ad  <- mean(b$SNP %chin% ad_nochr6_hdl$SNP)
  data.table(CHR = chr, piece = piece, n_ref = n_ref, cov_scz = cov_scz, cov_ad = cov_ad)
}

with_progress({
  p <- progressor(steps = nrow(NEWLOC_no6))
  cov_list_noChr6 <- lapply(seq_len(nrow(NEWLOC_no6)), function(i) {
    p(sprintf("coverage noChr6 chr %d piece %d", NEWLOC_no6$CHR[i], NEWLOC_no6$piece[i]))
    locus_coverage_noChr6(NEWLOC_no6$CHR[i], NEWLOC_no6$piece[i])
  })
})

cov_tab_noChr6 <- rbindlist(cov_list_noChr6)
cov_min_noChr6 <- 0.90
idx_all  <- seq_len(nrow(NEWLOC_no6))
idx_keep <- idx_all[(cov_tab_noChr6$cov_scz >= cov_min_noChr6) & (cov_tab_noChr6$cov_ad >= cov_min_noChr6)]
if (length(idx_keep) == 0L) idx_keep <- idx_all
message(sprintf("[noChr6] Loci kept by coverage ≥ %.0f%%: %d / %d",
                100*cov_min_noChr6, length(idx_keep), nrow(NEWLOC_no6)))

runner_noChr6 <- function(chr, piece) {
  HDL.L(
    gwas1 = sz_nochr6_hdl, gwas2 = ad_nochr6_hdl,
    Trait1name = "SCZ", Trait2name = "AD",
    LD.path = LD.path, bim.path = bim.path,
    chr = chr, piece = piece,
    eigen.cut = 0.99, N0 = 0, lim = exp(-18), output.file = ""
  )
}

with_progress({
  p2 <- progressor(steps = length(idx_keep))
  res_list_noChr6 <- lapply(idx_keep, function(i) {
    chr   <- NEWLOC_no6$CHR[i]
    piece <- NEWLOC_no6$piece[i]
    p2(sprintf("HDL-L noChr6 chr %d piece %d", chr, piece))
    try(runner_noChr6(chr, piece), silent = TRUE)
  })
})

ok_noChr6 <- vapply(res_list_noChr6, function(x) !(inherits(x, "try-error") || is.null(x)), logical(1))
dir.create("HDLL_out_noChr6", showWarnings = FALSE)
if (!any(ok_noChr6)) {
  fwrite(cov_tab_noChr6, "HDLL_out_noChr6/SCZ_AD.noChr6.HDLL.coverage_only.tsv", sep = "\t")
  stop("[noChr6] No loci returned valid HDL-L results. Wrote coverage table for debugging.")
}
res_local_noChr6 <- rbindlist(res_list_noChr6[ok_noChr6], fill = TRUE)

if (all(c("chr","piece") %in% names(res_local_noChr6))) {
  setnames(res_local_noChr6, c("chr","piece"), c("CHR","piece"))
}
res_join_noChr6 <- merge(res_local_noChr6, cov_tab_noChr6, by = c("CHR","piece"), all.x = TRUE)
if (all(c("CHR","piece") %in% names(res_join_noChr6))) {
  setnames(res_join_noChr6, c("CHR","piece"), c("chr","piece"))
}

fwrite(res_join_noChr6, "HDLL_out_noChr6/SCZ_AD.noChr6.HDLL.local_rg.with_coverage.tsv", sep = "\t")
saveRDS(res_join_noChr6, "HDLL_out_noChr6/SCZ_AD.noChr6.HDLL.local_rg.with_coverage.rds")
print(dim(res_join_noChr6))

valid_noChr6 <- res_join_noChr6[
  is.finite(Heritability_1) & is.finite(Heritability_2) & is.finite(Genetic_Covariance) &
    Heritability_1 > 0 & Heritability_2 > 0 &
    cov_scz >= cov_min_noChr6 & cov_ad >= cov_min_noChr6
]

if (nrow(valid_noChr6) > 1L) {
  h1   <- valid_noChr6$Heritability_1
  h2   <- valid_noChr6$Heritability_2
  gcov <- valid_noChr6$Genetic_Covariance
  
  rg_hat <- sum(gcov) / sqrt(sum(h1) * sum(h2))
  
  theta_i <- sapply(seq_len(nrow(valid_noChr6)), function(i) {
    h1_i   <- sum(h1) - h1[i]
    h2_i   <- sum(h2) - h2[i]
    gcov_i <- sum(gcov) - gcov[i]
    if (h1_i <= 0 || h2_i <= 0) return(NA_real_)
    gcov_i / sqrt(h1_i * h2_i)
  })
  theta_i <- theta_i[is.finite(theta_i)]
  n <- length(theta_i)
  theta_dot <- mean(theta_i)
  se_jk <- sqrt(((n - 1) / n) * sum((theta_i - theta_dot)^2))
  z <- rg_hat / se_jk
  p <- 2 * pnorm(-abs(z))
  ci_low  <- max(-1, rg_hat - 1.96 * se_jk)
  ci_high <- min(1, rg_hat + 1.96 * se_jk)
  
  global_tab_noChr6 <- data.table(
    rg = rg_hat, SE = se_jk, z = z, p = p,
    CI_low = ci_low, CI_high = ci_high,
    n_loci_used = nrow(valid_noChr6),
    cov_min_used = cov_min_noChr6
  )
  fwrite(global_tab_noChr6, "HDLL_out_noChr6/SCZ_AD.noChr6.HDLL.global_rg.tsv", sep = "\t")
  saveRDS(global_tab_noChr6, "HDLL_out_noChr6/SCZ_AD.noChr6.HDLL.global_rg.rds")
  print(global_tab_noChr6)
} else {
  message("[noChr6] Not enough valid loci (positive h2 & finite gcov under coverage rule) to compute a global rg.")
}
