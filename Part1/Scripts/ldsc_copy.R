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
remotes::install_github("josefin-werme/LAVA")

library(GenomicSEM)
packageVersion("GenomicSEM")

library(LAVA)
packageVersion("LAVA")

# HDL-L package - check for genetic correlation with this one vs LDSC
remotes::install_github("zhenin/HDL/HDL", dependencies = TRUE)

library(HDL)
library(dplyr)
library(data.table)
library(parallel)
packageVersion("HDL")

# You already did this earlier, but do it again to be explicit:
setwd("/Users/guillermocomesanacimadevila/Desktop/PhD/Part1/Data")

# 1) Run munge (as you have)
munged <- munge(
  files       = c("SZ/SZ.ldsc.sumstats", "AD/AD.ldsc.sumstats"),
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
# Run 4 -> Change INFO from 0.75 to 0.90 

# Run 5 -> rg = 0.08258438, SE = 0.04301183, z = 1.920, p = 0.05485
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

# rg: 0.1363 
# CI: (-0.9602, 1)
# P = 0.84

# WHOLE GENOME NOW 
# 1. Make sure == numeric
setDT(sz_gwas); setDT(ad_gwas)
sz_gwas[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]
ad_gwas[, `:=`(N = as.numeric(N), Z = as.numeric(Z))]

# 2. Run
runner <- function(chr, piece) try(
  HDL.L(
    gwas1 = sz_gwas, gwas2 = ad_gwas,
    Trait1name = "SCZ", Trait2name = "AD",
    LD.path = LD.path, bim.path = bim.path,
    chr = as.numeric(chr), piece = as.numeric(piece),
    eigen.cut = 0.99, N0 = 0, lim = exp(-18), output.file = ""
  ),
  silent = TRUE
)

# Slow as hell - Use * CPU cores, except 1 
cores <- max(1, detectCores() - 1)
res_list <- mcmapply(runner, NEWLOC$CHR, NEWLOC$piece, SIMPLIFY = FALSE, mc.cores = cores)
ok <- vapply(res_list, function(x) !(inherits(x,"try-error") || is.null(x)), logical(1))
res_local <- rbindlist(res_list[ok], fill = TRUE)

dir.create("HDLL_out", showWarnings = FALSE)
fwrite(res_local, "HDLL_out/SCZ_AD.HDLL.local_rg.csv", sep = ",")
saveRDS(res_local, "HDLL_out/SCZ_AD.HDLL.local_rg.rds")


setDT(res_local)
if ("Genetic.Correlation" %in% names(res_local)) res_local[, rg := Genetic.Correlation]
pcol <- if ("P.GC" %in% names(res_local)) "P.GC" else if ("P" %in% names(res_local)) "P" else NA_character_
if (!is.na(pcol)) {
  res_local[, p_rg := as.numeric(get(pcol))]
  res_local[, fdr := p.adjust(p_rg, method = "BH")]
  fwrite(res_local[order(fdr)], "HDLL_out/SCZ_AD.HDLL.local_rg.withFDR.tsv", sep = "\t")
}

print(dim(res_local))

# Weighted average total for HDL-L vs LDSC -> compare it 


