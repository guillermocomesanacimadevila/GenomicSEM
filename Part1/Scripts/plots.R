#!/usr/bin/env Rscript
# Change y/axis upper limit to 35


suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("qqman",       quietly = TRUE)) install.packages("qqman")
  library(data.table); library(qqman)
})

# -- helpers ---------------------------------------------------------------
with_alpha <- function(cols, a = 0.80) {
  vapply(cols, function(c) grDevices::adjustcolor(c, alpha.f = a), "", USE.NAMES = FALSE)
}

make_synth <- function(n_per_chr = 4000, seed = 1, n_hits = 60) {
  set.seed(seed)
  CHR  <- rep(1:22, each = n_per_chr)
  lens <- round(seq(60e6, 250e6, length.out = 22))
  BP   <- unlist(lapply(seq_along(lens), function(i) sort(sample.int(lens[i], n_per_chr))))
  P    <- runif(length(CHR))
  hit_idx <- sample.int(length(P), size = n_hits)
  P[hit_idx] <- runif(n_hits, min = 1e-12, max = 5e-8)
  data.frame(SNP = paste0("rs", seq_along(P)), CHR = CHR, BP = BP, P = P)
}

# -- plotting --------------------------------------------------------------
plot_two_manhattans <- function(
    d_top, d_bottom,
    out    = "~/Downloads/manhattan_ad_top_scz_bottom.nature.pdf",
    width  = 7.5,
    height = 5.8
) {
  # palette + thresholds
  blues   <- c("#0B3C5D", "#1F5A89", "#2F7FBF", "#66ACDC")
  blues_a <- with_alpha(blues, 0.80)
  gw_p  <- 5e-8;  sug_p <- 1e-5
  gw  <- -log10(gw_p);  sug <- -log10(sug_p)
  
  # fixed y-limit
  ymax <- 12
  
  # PDF device (no XQuartz needed)
  try(graphics.off(), silent = TRUE)
  out_path <- path.expand(out)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  pdf(out_path, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  
  # base style
  op <- par(
    mfrow = c(2, 1),
    mgp   = c(2.0, 0.6, 0),
    oma   = c(0, 0, 0, 0),    # remove outer margins
    mar   = c(1.0, 3.8, 1.2, 0.8), # will adjust per panel below
    xaxs  = "i",
    tcl   = -0.25,
    lwd   = 1.0,
    lend  = "butt",
    bty   = "l"
  )
  on.exit(par(op), add = TRUE)
  
  add_overlays <- function() {
    usr <- par("usr")
    ybreaks <- pretty(c(0, usr[4]), n = 6)
    abline(h = ybreaks, col = "#E0E0E0", lwd = 0.6)
    abline(h = gw,  col = "#7D3C98", lty = 2, lwd = 1.6)
    abline(h = sug, col = "#BFC5CA", lty = 2, lwd = 1.2)
  }
  
  # -------- TOP panel (A): NO CHR LABELS --------
  par(mar = c(0.8, 3.8, 1.0, 0.8))  # tighten spacing below
  qqman::manhattan(
    d_top,
    col            = blues_a,
    ylim           = c(0, ymax),
    genomewideline = FALSE,
    suggestiveline = FALSE,
    cex            = 0.5,
    cex.axis       = 0.9,
    cex.lab        = 1.0,
    ylab           = expression(-log[10](P)),
    xlab           = "",
    chrlabs        = rep("", length(unique(d_top$CHR))),  # hide chr labels
    main           = ""
  )
  add_overlays()
  usr <- par("usr")
  text(x = usr[1] - 0.3 * diff(usr[1:2]), y = ymax * 0.9,
       labels = "A", font = 2, cex = 1.2, xpd = NA)
  
  # -------- BOTTOM panel (B): KEEP CHR LABELS --------
  par(mar = c(3.0, 3.8, 0.8, 0.8))  # tighten gap above
  qqman::manhattan(
    d_bottom,
    col            = blues_a,
    ylim           = c(0, ymax),
    genomewideline = FALSE,
    suggestiveline = FALSE,
    cex            = 0.5,
    cex.axis       = 0.9,
    cex.lab        = 1.0,
    ylab           = expression(-log[10](P)),
    main           = ""
  )
  add_overlays()
  usr <- par("usr")
  text(x = usr[1] - 0.3 * diff(usr[1:2]), y = ymax * 0.9,
       labels = "B", font = 2, cex = 1.2, xpd = NA)
  
  invisible(normalizePath(out_path))
}

# -- demo with synthetic data ---------------------------------------------
d_ad <- make_synth(n_per_chr = 4000, seed = 101, n_hits = 80)  # TOP
d_sz <- make_synth(n_per_chr = 4000, seed = 202, n_hits = 60)  # BOTTOM

outfile <- plot_two_manhattans(
  d_ad, d_sz,
  out = "~/Downloads/manhattan_ad_top_scz_bottom.nature.pdf"
)

cat("✔ Saved:", outfile, "\n")
