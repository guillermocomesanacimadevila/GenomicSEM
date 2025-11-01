#!/usr/bin/env Rscript
# Running LDSC (AD -> SCZ) 

suppressPackageStartupMessages({library(GenomicSEM);library(data.table)})

proj <- "/Users/c24102394/Desktop/PhD/AD_SCZ_AGE"
data_dir <- file.path(proj, "Data")
setwd(data_dir)

hm3_path <- "/Users/c24102394/ref/ldsc/w_hm3.snplist"
ld_path <- "/Users/c24102394/ref/ldsc/eur_w_ld_chr"
wld_path <- "/Users/c24102394/ref/ldsc/weights_hm3_no_hla"
in_files <- c(
  file.path(data_dir, "AD/post-qc/Kunkle_2019_IGAP.ldsc_ready_neff.tsv"),
  file.path(data_dir, "SCZ/post-qc/PGC3_SCZ_wave3.ldsc_ready_neff.tsv")
)

# First AD-SCZ
trait_names <- c("AD", "SCZ")
N_vec <- c(AD=88394, SCZ=157013)

munge(
  files = in_files,
  maf.filter = 0.01,
  info.filter = 0.9,
  N = as.numeric(N_vec[trait_names]),
  hm3 = hm3_path,
  trait.names = trait_names
)

out_files <- file.path(getwd(),paste0(trait_names,".sumstats.gz"))
stopifnot(all(file.exists(out_files)))

ad_cases <- 35274; ad_controls <- 59163
scz_cases <- 67390; scz_controls <- 94015
sample_prev <- c(
  AD = ad_cases / (ad_cases+ad_controls),
  SCZ = scz_cases / (scz_cases+scz_controls)
)
population_prev <- c(AD=0.07, SCZ=0.01)
ldsc_out <- ldsc(
  traits=out_files,
  sample.prev = as.numeric(sample_prev[trait_names]),
  population.prev = as.numeric(population_prev[trait_names]),
  ld = ld_path,
  wld = wld_path,
  ldsc.log = TRUE
)

uni_AD <- ldsc(
  traits = out_files[1],
  ld = ld_path,
  wld = wld_path,
  population.prev = as.numeric(population_prev[1]),
  sample.prev = as.numeric(sample_prev[1])
)

uni_SCZ <- ldsc(
  traits = out_files[2],
  ld = ld_path,
  wld = wld_path,
  population.prev = as.numeric(population_prev[2]),
  sample.prev = as.numeric(sample_prev[2])
)

get_mean <- function(obj) {
  df <- as.data.frame(obj$ldscoutput)
  nms <- names(df)
  k <- grep("Mean.*Chi", nms, ignore.case=TRUE)
  if (length(k)) as.numeric(df[[k[1]]]) else NA_real_
}

# ldsc_out
# uni_AD
# uni_BIP

qc_out <- rbindlist(list(
  data.table(
    Trait ="AD",
    h2 = as.numeric(uni_AD$S[1,1]),
    h2_se = sqrt(as.numeric(uni_AD$V[1,1])),
    mean_chisq = get_mean(uni_AD),
    Intercept = as.numeric(uni_AD$I[1,1])
  ),
  data.table(
    Trait = "SCZ",
    h2 = as.numeric(uni_SCZ$S[1,1]),
    h2_se = sqrt(as.numeric(uni_SCZ$V[1,1])),
    mean_chisq = get_mean(uni_SCZ),
    Intercept = as.numeric(uni_SCZ$I[1,1])
  )
))

# qc_out
qc_out[, h2_z := h2/h2_se]
qc_out[, pass_h2z := h2_z > 2]
qc_out[, pass_mean := ifelse(is.na(mean_chisq), NA, mean_chisq>1.02)]
qc_out[, pass_int := Intercept>=0.9 & Intercept<=1.1]
qc_out[, pass_all := ifelse(is.na(pass_mean), pass_h2z & pass_int, pass_h2z & pass_mean & pass_int)]
print(qc_out[,.(Trait,h2,h2_se,h2_z,mean_chisq,Intercept,pass_all)])

outdir <- file.path(proj,"outputs/ldsc/ad_scz"); dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
fwrite(qc_out,file.path(outdir,"ldsc_trait_qc_AD_SCZ.tsv"),sep="\t")
write.csv(ldsc_out$S,file.path(outdir,"ldsc_S.csv"),row.names=FALSE)
write.csv(ldsc_out$V,file.path(outdir,"ldsc_V.csv"),row.names=FALSE)
write.csv(ldsc_out$N,file.path(outdir,"ldsc_N.csv"),row.names=FALSE)
write.csv(data.frame(m=ldsc_out$m),file.path(outdir,"ldsc_m.csv"),row.names=FALSE)

ldsc_rg_summary <- function(ld) {
  S <- ld$S
  V <- ld$V
  rg <- S[1, 2] / sqrt(S[1, 1] * S[2, 2])
  g <- c(
    -0.5 * rg / S[1, 1],
    1 / sqrt(S[1, 1] * S[2, 2]),
    -0.5 * rg / S[2, 2]
  )
  se <- sqrt(as.numeric(t(g) %*% V %*% g))
  z <- rg / se
  p <- 2 * pnorm(-abs(z))
  zc <- qnorm(0.975)
  ci <- c(
    max(-1, rg - zc * se),
    min(1, rg + zc * se)
  )
  data.frame(
    rg = rg,
    SE = se,
    z = z,
    p = p,
    CI_low = ci[1],
    CI_high = ci[2]
  )
}

print(ldsc_rg_summary(ldsc_out))

# Build sample overlap correlation matrix 
I_mat <- ldsc_out$I
rownames(I_mat) <- colnames(I_mat) <- c("AD","SCZ")
overlap_corr <- I_mat / sqrt(outer(diag(I_mat), diag(I_mat), "*"))
diag(overlap_corr) <- 1
write.csv(overlap_corr, file.path(outdir, "overlap_corr_for_LAVA_AD_SCZ.csv"), quote=FALSE)
print(overlap_corr)