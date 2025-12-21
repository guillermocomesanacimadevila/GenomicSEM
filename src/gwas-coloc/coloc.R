#!/usr/bin/env Rscript
# FOR NEXTFLOW PROCESSES
# mkdir -p ${workflow.launchDir}/outputs/coloc/AD_SCZ
# mkdir -p ${workflow.launchDir}/outputs/coloc/AD_LONG
# mkdir -p ${workflow.launchDir}/outputs/coloc/SCZ_LONG

suppressPackageStartupMessages({
  library(coloc)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop(
    "Usage: Rscript coloc.R ",
    "<prefix> <trait1> <trait2> <loci_dir> <out_dir> <type1> <type2> [s1] [s2]\n",
    "prefix = AD_SCZ or AD_LONG_SCZ_LONG\n",
    "trait1 = AD / SCZ / LONG\n",
    "trait2 = AD / SCZ / LONG\n",
    "loci_dir = dir with <prefix>_locus_<id>_<trait>.tsv\n",
    "out_dir = dir for results\n",
    "type1 = cc or quant\n",
    "type2 = cc or quant\n",
    "s1,s2 = case fractions for cc traits (ignored for quant)\n"
  )
}

prefix <- args[1]
trait1 <- args[2]
trait2 <- args[3]
loci_dir <- args[4]
out_dir <- args[5]
type1 <- args[6]
type2 <- args[7]

if (length(args) >= 8) {
  s1 <- as.numeric(args[8])
} else {
  s1 <- 0.5
}
if (length(args) >= 9) {
  s2 <- as.numeric(args[9])
} else {
  s2 <- 0.5
}

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
out_file <- file.path(out_dir, paste(prefix, trait1, trait2, "coloc.tsv", sep="_"))

pattern1 <- paste0("(^", prefix, "_)?locus_[0-9]+_", trait1, "\\.tsv$")
pattern2 <- paste0("(^", prefix, "_)?locus_[0-9]+_", trait2, "\\.tsv$")
all_files <- list.files(loci_dir, pattern="locus_.*\\.tsv$", full.names=TRUE)
base_names <- basename(all_files)
files1 <- all_files[grepl(pattern1, base_names)]
files2 <- all_files[grepl(pattern2, base_names)]
get_id <- function(f) {
  bn <- basename(f)
  parts <- strsplit(bn, "_")[[1]]
  if (startsWith(bn, "locus_")) {
    return(as.integer(parts[2]))
  }
  as.integer(parts[3])
}

if (length(files1)==0 || length(files2)==0) {
  stop("No locus files found for this prefix/trait combination: ",
       "prefix=", prefix, " trait1=", trait1, " trait2=", trait2)
}

df1 <- data.frame(id=sapply(files1, get_id), f1=files1, stringsAsFactors=FALSE)
df2 <- data.frame(id=sapply(files2, get_id), f2=files2, stringsAsFactors=FALSE)

pairs <- merge(df1, df2, by="id")
if (nrow(pairs)==0) {
  stop("No matching loci between ", trait1, " and ", trait2, ".")
}

cat("Running coloc for", nrow(pairs), "loci (", trait1, " vs ", trait2, ")\n", sep="")

results <- list()

for (i in seq_len(nrow(pairs))) {
  locus_id <- pairs$id[i]
  f1 <- pairs$f1[i]
  f2 <- pairs$f2[i]
  d1 <- read.table(f1, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  d2 <- read.table(f2, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  m <- merge(d1, d2, by="SNP", suffixes=c(".1", ".2"))
  same <- m$A1.1 == m$A1.2 & m$A2.1 == m$A2.2
  flip <- m$A1.1 == m$A2.2 & m$A2.1 == m$A1.2
  m$BETA.2[flip] <- -m$BETA.2[flip]
  m <- m[same | flip, ]
  lead_snp <- NA_character_
  if (nrow(m) > 0) {
    if ("P.1" %in% names(m)) {
      idx <- which.min(m$P.1)
      lead_snp <- m$SNP[idx]
    } else if ("P.2" %in% names(m)) {
      idx <- which.min(m$P.2)
      lead_snp <- m$SNP[idx]
    }
  }
  
  get_N1 <- function(mm) {
    if ("N.1" %in% names(mm)) return(mm$N.1)
    if ("N_eff.1" %in% names(mm)) return(mm$N_eff.1)
    if ("Neff.1" %in% names(mm)) return(mm$Neff.1)
    if ("NEFF.1" %in% names(mm)) return(mm$NEFF.1)
    return(rep(NA_real_, nrow(mm)))
  }
  
  get_N2 <- function(mm) {
    if ("N.2" %in% names(mm)) return(mm$N.2)
    if ("N_eff.2" %in% names(mm)) return(mm$N_eff.2)
    if ("Neff.2" %in% names(mm)) return(mm$Neff.2)
    if ("NEFF.2" %in% names(mm)) return(mm$NEFF.2)
    return(rep(NA_real_, nrow(mm)))
  }
  
  if (nrow(m) < 10) {
    res_row <- data.frame(
      nsnps = NA,
      PP.H0.abf = NA,
      PP.H1.abf = NA,
      PP.H2.abf = NA,
      PP.H3.abf = NA,
      PP.H4.abf = NA,
      locus_id = locus_id,
      lead_snp = lead_snp,
      n_snps = nrow(m),
      stringsAsFactors=FALSE
    )
  } else {
    ds1 <- list(
      snp = m$SNP,
      beta = m$BETA.1,
      varbeta = m$SE.1^2,
      N = get_N1(m),
      type = type1
    )
    if (type1 == "cc") {
      ds1$s <- s1
    }
    if (type1 == "quant") {
      ds1$sdY <- 1
    }
    
    ds2 <- list(
      snp = m$SNP,
      beta = m$BETA.2,
      varbeta = m$SE.2^2,
      N = get_N2(m),
      type = type2
    )
    if (type2 == "cc") {
      ds2$s <- s2
    }
    if (type2 == "quant") {
      ds2$sdY <- 1
    }
    
    co <- coloc.abf(ds1, ds2)
    s  <- as.data.frame(t(co$summary), stringsAsFactors=FALSE)
    s$locus_id <- locus_id
    s$lead_snp <- lead_snp
    s$n_snps   <- nrow(m)
    res_row <- s
  }
  
  results[[length(results)+1]] <- res_row
  pct <- round(100*i/nrow(pairs))
  width  <- 40
  filled <- round(width*i/nrow(pairs))
  bar <- paste0(
    paste(rep("#", filled), collapse=""),
    paste(rep("-", width-filled), collapse="")
  )
  cat(sprintf("\r[%s] %d/%d (%d%%)", bar, i, nrow(pairs), pct))
  flush.console()
}

cat("\nDone.\n")

all_cols <- unique(unlist(lapply(results, names)))
results2 <- lapply(results, function(x) {
  missing <- setdiff(all_cols, names(x))
  if (length(missing) > 0) {
    for (nm in missing) {
      x[[nm]] <- NA
    }
  }
  x <- x[all_cols]
  x
})

final <- do.call(rbind, results2)
if ("locus_id" %in% names(final)) {
  final <- final[order(final$locus_id), ]
}

write.table(final, out_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote:", out_file, "\n")