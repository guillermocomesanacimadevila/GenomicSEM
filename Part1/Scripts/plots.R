#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman")
  library(data.table); library(qqman)
})

# ---- synthetic data generator ----
make_synth <- function(n_per_chr = 4000, seed = 1, n_hits = 60) {
  set.seed(seed)
  CHR <- rep(1:22, each = n_per_chr)
  lens <- round(seq(60e6, 250e6, length.out = 22))
  BP <- unlist(lapply(seq_along(lens), function(i) sort(sample.int(lens[i], n_per_chr))))
  P <- runif(length(CHR))
  hit_idx <- sample.int(length(P), size = n_hits)
  P[hit_idx] <- runif(n_hits, min = 1e-12, max = 5e-8)
  data.frame(SNP = paste0("rs", seq_along(P)), CHR = CHR, BP = BP, P = P)
}

# ---- plotting function ----
plot_two_manhattans <- function(d_top, d_bottom,
                                out = "~/Downloads/manhattan_ad_top_scz_bottom.synthetic.pdf",
                                width = 10, height = 8) {
  
  # blue palette
  blues <- c("#0B3C5D", "#1F5A89", "#2F7FBF", "#66ACDC")
  gw  <- -log10(5e-8)
  sug <- -log10(1e-5)
  
  # safely open PDF device
  try(graphics.off(), silent = TRUE)
  out_path <- path.expand(out)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(capabilities("cairo"))) {
    cairo_pdf(out_path, width = width, height = height)
  } else {
    pdf(out_path, width = width, height = height)
  }
  on.exit(dev.off(), add = TRUE)
  
  op <- par(mfrow = c(2,1), mgp = c(2.2, 0.6, 0))
  on.exit(par(op), add = TRUE)
  
  # --- TOP panel (A)
  par(mar = c(1.8, 4.2, 2.2, 0.8))
  qqman::manhattan(
    d_top,
    main = "",
    col  = blues,
    genomewideline = FALSE,
    suggestiveline = FALSE,
    cex = 0.45, pch = 20,
    cex.axis = 0.95, cex.lab = 1.0,
    ylab = expression(-log[10](P))
  )
  abline(h = gw,  col = "#8E44AD", lty = 2, lwd = 2)
  abline(h = sug, col = "#BDBDBD", lty = 2, lwd = 1.5)
  mtext("A", side = 3, adj = 0, line = 0.2, font = 2, cex = 1.1)
  
  # --- BOTTOM panel (B)
  par(mar = c(4.6, 4.2, 2.2, 0.8))
  qqman::manhattan(
    d_bottom,
    main = "",
    col  = blues,
    genomewideline = FALSE,
    suggestiveline = FALSE,
    cex = 0.45, pch = 20,
    cex.axis = 0.95, cex.lab = 1.0,
    ylab = expression(-log[10](P))
  )
  abline(h = gw,  col = "#8E44AD", lty = 2, lwd = 2)
  abline(h = sug, col = "#BDBDBD", lty = 2, lwd = 1.5)
  mtext("B", side = 3, adj = 0, line = 0.2, font = 2, cex = 1.1)
  
  invisible(normalizePath(out_path))
}

# ---- run with synthetic data ----
d_ad <- make_synth(n_per_chr = 4000, seed = 101, n_hits = 80)  # AD (top)
d_sz <- make_synth(n_per_chr = 4000, seed = 202, n_hits = 60)  # SCZ (bottom)

outfile <- plot_two_manhattans(d_ad, d_sz,
                               out = "~/Downloads/manhattan_ad_top_scz_bottom.synthetic.pdf")

cat("Saved:", outfile, "\n")
