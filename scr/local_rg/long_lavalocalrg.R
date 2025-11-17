#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(LAVA); library(progress); library(data.table)})

proj <- "/Users/c24102394/Desktop/PhD/AD_SCZ_AGE"
data_dir <- file.path(proj,"Data")
REF_PREFIX <- "/Users/c24102394/Desktop/PhD/GenomicSEM/Part1/Checks/Checks/lava_ref/lava-ukb-v1.1"
LOCI_FILE <- "/Users/c24102394/Desktop/PhD/AD_SZ_genes/LAVA_ref/hdll_blocks.coords.loci"
OUT_DIR <- file.path(proj, "outputs", "lava/ad_scz_age"); dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
setwd(data_dir)

in_files <- c(
  file.path(data_dir,"AD/post-lava/AD_lava_ready.tsv"),
  file.path(data_dir,"AGE/post-lava/AGE_lava_ready.tsv"),
  file.path(data_dir,"SCZ/post-lava/SCZ_lava_ready.tsv")
)

PHENOS <- c("AD","AGE","SCZ")
info <- data.frame(
  phenotype = PHENOS,
  cases     = c(35274, 354854, 67390),
  controls  = c(59163, 354855, 94015),
  prevalence= c(0.07, 0.10, 0.01),
  filename  = in_files,
  stringsAsFactors = FALSE
)

INFO_TSV <- file.path(proj, "InputFiles/input.info_ad_bip_age.txt")
dir.create(dirname(INFO_TSV), recursive=TRUE, showWarnings=FALSE)
fwrite(info, INFO_TSV, sep="\t", quote=FALSE, na="NA")

loci <- read.loci(LOCI_FILE)
p_gate <- 0.05 / nrow(loci)

is_mhc <- function(chr, start, stop){
  chr6  <- as.character(chr) %in% c("6","chr6")
  inwin <- !(stop < 25000000 | start > 34000000)
  isTRUE(chr6 & inwin)
}

ov_ab <- as.matrix(read.csv(file.path(proj, "outputs/ldsc/ad_age/overlap_corr_for_LAVA_AD_AGE.csv"), row.names = 1, check.names = FALSE))
ov_as <- as.matrix(read.csv(file.path(proj, "outputs/ldsc/ad_scz/overlap_corr_for_LAVA_AD_SCZ.csv"), row.names = 1, check.names = FALSE))
ov_sb <- as.matrix(read.csv(file.path(proj, "outputs/ldsc/scz_age/overlap_corr_for_LAVA_SCZ_AGE.csv"), row.names = 1, check.names = FALSE))
OV <- diag(3)
dimnames(OV) <- list(PHENOS, PHENOS)
OV["AD","AGE"] <- OV["AGE","AD"] <- ov_ab["AD","AGE"]
OV["AD","SCZ"] <- OV["SCZ","AD"] <- ov_as["AD","SCZ"]
OV["SCZ","AGE"] <- OV["AGE","SCZ"] <- ov_sb["SCZ","AGE"]
overlap_file <- file.path(proj, "outputs/ldsc/overlap_corr_for_LAVA_AD_AGE_SCZ.csv")
write.table(OV, overlap_file, sep = "\t", quote = FALSE, col.names = NA)

inp <- process.input(
  input.info.file     = INFO_TSV,
  sample.overlap.file = overlap_file,
  ref.prefix          = REF_PREFIX,
  phenos              = PHENOS
)

pairs <- list(c("AD","AGE"), c("AD","SCZ"), c("SCZ","AGE"))
pb <- progress_bar$new(total=nrow(loci), width=60)
rows_bivar <- vector("list", nrow(loci)*length(pairs))
rows_partial <- vector("list", nrow(loci)*length(pairs))
ri <- 0L; rj <- 0L

for (i in seq_len(nrow(loci))) {
  pb$tick()
  loc <- process.locus(loci[i,], inp)
  ph <- intersect(loc$phenos, PHENOS)
  if (length(ph)==0) next
  uni <- run.univ(loc, phenos=ph)
  h2 <- setNames(rep(NA_real_, length(PHENOS)), PHENOS)
  pu <- setNames(rep(NA_real_, length(PHENOS)), PHENOS)
  if (!is.null(uni) && nrow(uni)>0) {
    for (phn in ph) {
      h2[phn] <- uni$h2.obs[uni$phen==phn]
      pu[phn] <- uni$p[uni$phen==phn]
    }
  }
  for (pp in pairs) {
    t1 <- pp[1]; t2 <- pp[2]
    if (!(t1 %in% ph && t2 %in% ph)) next
    gate1 <- is.finite(h2[t1]) && h2[t1]>0 && is.finite(pu[t1]) && pu[t1] < p_gate
    gate2 <- is.finite(h2[t2]) && h2[t2]>0 && is.finite(pu[t2]) && pu[t2] < p_gate
    if (gate1 && gate2) {
      bv <- try(run.bivar(loc, phenos=c(t1,t2), param.lim=2), silent=TRUE)
      if (!inherits(bv,"try-error") && !is.null(bv) && nrow(bv)>0) {
        lc <- if (is.finite(bv$rho) && is.finite(h2[t1]) && is.finite(h2[t2]) && h2[t1]>=0 && h2[t2]>=0) bv$rho*sqrt(h2[t1]*h2[t2]) else NA_real_
        ri <- ri + 1L
        rows_bivar[[ri]] <- data.frame(
          locus=loc$id, chr=loc$chr, start=loc$start, stop=loc$stop,
          pair=paste(t1,t2,sep="_"),
          phen1=t1, phen2=t2,
          n_snps=if (!is.null(loc$n.snps)) as.integer(loc$n.snps) else NA_integer_,
          h2_1=h2[t1], p_1=pu[t1], h2_2=h2[t2], p_2=pu[t2],
          rho=bv$rho, rho.lower=bv$rho.lower, rho.upper=bv$rho.upper, r2=bv$r2, p=bv$p,
          local_cov=lc, stringsAsFactors=FALSE
        )
        third <- setdiff(PHENOS, c(t1,t2))
        if (length(third)==1 && third %in% ph) {
          gate3 <- is.finite(h2[third]) && h2[third]>0 && is.finite(pu[third]) && pu[third] < p_gate
          if (gate3) {
            par <- try(run.pcor(loc, phenos=ph, target=c(t1,t2), param.lim=2), silent=TRUE)
            if (!inherits(par,"try-error") && !is.null(par) && nrow(par)>0) {
              pr <- par[par$phen1==t1 & par$phen2==t2, , drop=FALSE]
              if (nrow(pr)==1) {
                rj <- rj + 1L
                rows_partial[[rj]] <- data.frame(
                  locus=loc$id, chr=loc$chr, start=loc$start, stop=loc$stop,
                  pair=paste(t1,t2,sep="_"),
                  phen1=t1, phen2=t2,
                  pcor=pr$pcor, p.partial=pr$p,
                  stringsAsFactors=FALSE
                )
              }
            }
          }
        }
      }
    }
  }
}

bivar <- rbindlist(rows_bivar[seq_len(ri)], use.names=TRUE, fill=TRUE)
partial <- rbindlist(rows_partial[seq_len(rj)], use.names=TRUE, fill=TRUE)

if (nrow(bivar)>0) {
  bivar[, mhc := is_mhc(as.integer(chr), start, stop)]
  idx <- which(is.finite(bivar$p))
  K_global <- length(idx)
  if (K_global>0) {
    alpha <- 0.05 / K_global
    bivar$p_bonf_paper <- NA_real_
    bivar$sig_paper <- NA
    bivar$q_fdr <- NA_real_
    bivar$p_bonf_paper[idx] <- bivar$p[idx] * K_global
    bivar$sig_paper[idx] <- bivar$p[idx] < alpha
    bivar$q_fdr[idx] <- p.adjust(bivar$p[idx], method="fdr")
  }
}

if (nrow(partial)>0 && "sig_paper" %in% names(bivar)) {
  key <- bivar[sig_paper==TRUE, .(locus, pair)]
  setkey(key, locus, pair); setkey(partial, locus, pair)
  partial <- partial[key, nomatch=0]
}

if (nrow(partial)>0) {
  for (pp in unique(partial$pair)) {
    idx <- which(partial$pair==pp & is.finite(partial$p.partial))
    Kp <- length(idx)
    if (Kp>0) {
      alpha_p <- 0.05 / Kp
      partial$p_bonf_paper.partial[idx] <- partial$p.partial[idx] * Kp
      partial$sig_paper.partial[idx] <- partial$p.partial[idx] < alpha_p
      partial$q_fdr.partial[idx] <- p.adjust(partial$p.partial[idx], method="fdr")
    }
  }
}

fwrite(bivar, file.path(OUT_DIR, "LAVA_local_rg_bivariate.tsv"), sep="\t")
fwrite(partial, file.path(OUT_DIR, "LAVA_local_rg_partial.tsv"), sep="\t")

bivar_compact <- bivar[, .(locus,chr,start,stop,mhc,pair,phen1,phen2,n_snps,h2_1,p_1,h2_2,p_2,rho,rho.lower,rho.upper,r2,p,p_bonf_paper,sig_paper,q_fdr,local_cov)]
bivar_compact <- if (nrow(bivar) > 0) {
  bivar[, .(locus,chr,start,stop,mhc,pair,phen1,phen2,n_snps,h2_1,p_1,h2_2,p_2,
            rho,rho.lower,rho.upper,r2,p,p_bonf_paper,sig_paper,q_fdr,local_cov)]
} else {
  data.table(locus=character(), chr=integer(), start=integer(), stop=integer(), mhc=logical(),
             pair=character(), phen1=character(), phen2=character(), n_snps=integer(),
             h2_1=numeric(), p_1=numeric(), h2_2=numeric(), p_2=numeric(),
             rho=numeric(), rho.lower=numeric(), rho.upper=numeric(), r2=numeric(), p=numeric(),
             p_bonf_paper=numeric(), sig_paper=logical(), q_fdr=numeric(), local_cov=numeric())
}

partial_compact <- if (nrow(partial) > 0) {
  partial[, .(locus,chr,start,stop,pair,phen1,phen2,pcor,p.partial,
              p_bonf_paper.partial,sig_paper.partial,q_fdr.partial)]
} else {
  data.table(locus=character(), chr=integer(), start=integer(), stop=integer(),
             pair=character(), phen1=character(), phen2=character(),
             pcor=numeric(), p.partial=numeric(),
             p_bonf_paper.partial=numeric(), sig_paper.partial=logical(), q_fdr.partial=numeric())
}

fwrite(bivar_compact, file.path(OUT_DIR, "LAVA_local_rg_bivariate.compact.csv"))
fwrite(partial_compact, file.path(OUT_DIR, "LAVA_local_rg_partial.compact.csv"))

pairs_vec <- if (nrow(bivar) > 0) unique(bivar$pair) else character(0)
if (length(pairs_vec) > 0) {
  for (pp in pairs_vec) {
    sub <- bivar[pair == pp]
    if (nrow(sub) == 0) next
    cov_all <- ifelse(is.finite(sub$local_cov) & !sub$mhc, sub$local_cov, 0)
    h1_all  <- pmax(ifelse(is.finite(sub$h2_1) & !sub$mhc, sub$h2_1, 0), 0)
    h2_all  <- pmax(ifelse(is.finite(sub$h2_2) & !sub$mhc, sub$h2_2, 0), 0)
    cov_sum <- sum(cov_all); h1_sum <- sum(h1_all); h2_sum <- sum(h2_all)
    if (h1_sum == 0 || h2_sum == 0) next
    rg_all  <- cov_sum / sqrt(h1_sum * h2_sum)
    M <- nrow(sub)
    if (M == 0) next
    rg_i <- numeric(M)
    for (k in seq_len(M)) {
      cs <- cov_sum - cov_all[k]
      h1s <- h1_sum - h1_all[k]
      h2s <- h2_sum - h2_all[k]
      if (h1s <= 0 || h2s <= 0) { rg_i[k] <- NA_real_; next }
      rg_i[k] <- cs / sqrt(h1s * h2s)
    }
    rg_i <- rg_i[is.finite(rg_i)]
    if (length(rg_i) == 0) next
    mu <- mean(rg_i)
    varJ <- (length(rg_i)-1)/length(rg_i) * sum((rg_i - mu)^2)
    seJ <- sqrt(varJ)
    if (!is.finite(seJ) || seJ == 0) next
    zJ <- rg_all / seJ
    pJ <- 2*pnorm(-abs(zJ))
    cat(sprintf("%s aggregate rg (exclude MHC): rg=%.4f SE=%.4f z=%.2f p=%.3g\n", pp, rg_all, seJ, zJ, pJ))
  }
}

cat("Saved:\n")
cat(file.path(OUT_DIR, "LAVA_local_rg_bivariate.tsv"), "\n")
cat(file.path(OUT_DIR, "LAVA_local_rg_partial.tsv"), "\n")
cat(file.path(OUT_DIR, "LAVA_local_rg_bivariate.compact.csv"), "\n")
cat(file.path(OUT_DIR, "LAVA_local_rg_partial.compact.csv"), "\n")