#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(gprofiler2)
  library(data.table)
})

# args - gene set 
# outdir
# pheno1_prefix
# pheno2_prefix

args <- commandArgs(trailingOnly = FALSE)
if (length(args) < 3) stop("Usage: conjFDR.R <gene_list> <out_prefix1> <out_prefix2> <out_dir>")
gene_list <- args[1]
out_prefix1 <- args[2]
out_prefix2 <- args[3]
out_dir <- args[4]

genes <- c("ADAM10", "MINDY2")

gostres <- gost(
  query = genes2,
  organism = "hsapiens",      
  user_threshold = 0.05,     
  correction_method = "g_SCS" 
)

res <- gostres$result
res