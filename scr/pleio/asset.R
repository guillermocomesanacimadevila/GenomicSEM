#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ASSET)
  library(progress)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript other_Asset.R <input_tsv> <out_prefix>")

infile <- args[1]
out_prefix <- args[2]

OUT_DIR <- "../../outputs/ASSET"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Reading input file:", infile, "\n")
dat <- fread(infile)
cat("Initial number of rows:", nrow(dat), "\n")

dat <- dat[
  !is.na(BETA_AD) & !is.na(SE_AD) &
    !is.na(BETA_SCZ) & !is.na(SE_SCZ) &
    !is.na(BETA_AGE) & !is.na(SE_AGE)
]
cat("Rows after filtering:", nrow(dat), "\n")

if (nrow(dat) == 0) {
  stop("No variants left after filtering")
}

traits.lab <- c("AD", "SCZ", "AGE")
ncase <- c(35274, 50965, 354854)
ncntl <- c(59163, 68049, 354855)
cor.mat <- diag(3)

n <- nrow(dat)
chunk_size <- 50000L
n_chunks <- ceiling(n / chunk_size)

cat("Total SNPs:", n, "\n")
cat("Chunk size:", chunk_size, "\n")
cat("Number of chunks:", n_chunks, "\n")

pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)

res_list <- vector("list", n_chunks)

for (i in seq_len(n_chunks)) {
  start <- (i - 1L) * chunk_size + 1L
  end <- min(i * chunk_size, n)
  dat_i <- dat[start:end]
  
  snp.vars <- dat_i$SNP
  beta.hat <- as.matrix(dat_i[, .(BETA_AD, BETA_SCZ, BETA_AGE)])
  sigma.hat <- as.matrix(dat_i[, .(SE_AD, SE_SCZ, SE_AGE)])
  
  res_i <- h.traits(
    snp.vars = snp.vars,
    traits.lab = traits.lab,
    beta.hat = beta.hat,
    sigma.hat = sigma.hat,
    ncase = ncase,
    ncntl = ncntl,
    cor = cor.mat,
    side = 2,
    search = 2,
    meta = TRUE
  )
  
  hs_i <- h.summary(res_i)
  
  meta_tab <- as.data.table(hs_i[[1]])
  setnames(meta_tab, "Pvalue", "P_ASSET")
  
  sub_tab <- NULL
  if ("Subset.2sided" %in% names(hs_i)) {
    sub_tab <- as.data.table(hs_i[["Subset.2sided"]])
    if ("Pheno" %in% names(sub_tab)) {
      setnames(sub_tab, "Pheno", "BEST_SUBSET")
    }
  }
  
  out_i <- merge(
    dat_i[, .(
      SNP, A1, A2,
      BETA_AD, SE_AD, P_AD,
      BETA_SCZ, SE_SCZ, P_SCZ,
      BETA_AGE, SE_AGE, P_AGE
    )],
    meta_tab[, .(SNP, P_ASSET)],
    by = "SNP",
    all.x = TRUE
  )
  
  if (!is.null(sub_tab) && "BEST_SUBSET" %in% names(sub_tab)) {
    out_i <- merge(
      out_i,
      sub_tab[, .(SNP, BEST_SUBSET)],
      by = "SNP",
      all.x = TRUE
    )
  } else {
    out_i[, BEST_SUBSET := NA_character_]
  }
  
  res_list[[i]] <- out_i
  
  setTxtProgressBar(pb, i)
}

close(pb)

out <- rbindlist(res_list, use.names = TRUE)

outfile1 <- file.path(OUT_DIR, paste0(out_prefix, "_ASSET_results.tsv"))
fwrite(out, outfile1, sep = "\t")
cat("\nSaved full results to:", outfile1, "\n")

pleio <- out[
  !is.na(P_ASSET) &
    P_ASSET < 5e-8 &
    P_AD < 0.05 & P_SCZ < 0.05 & P_AGE < 0.05 &
    !is.na(BEST_SUBSET) &
    grepl("AD", BEST_SUBSET) &
    grepl("SCZ", BEST_SUBSET) &
    grepl("AGE", BEST_SUBSET)
]

outfile2 <- file.path(OUT_DIR, paste0(out_prefix, "_ASSET_pleiotropic_hits.tsv"))
fwrite(pleio, outfile2, sep = "\t")
cat("Saved pleiotropic hits to:", outfile2, "\n")

cat("All done.\n")
