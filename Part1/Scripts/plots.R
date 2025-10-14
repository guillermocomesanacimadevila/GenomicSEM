#!/usr/bin/env Rscript
# AD (top) + SCZ (bottom) Manhattan figure
# - y-axis top limit = 35
# - no chromosome labels on TOP panel
# - tighter spacing between panels
# - panel labels "A" and "B" near lower-left

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("qqman",       quietly = TRUE)) install.packages("qqman")
  library(data.table); library(qqman)
})

# --- paths (edit if your folders differ) ----------------------------------
# AD (Kunkle) file:
AD_FILE <- "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/AD/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.prepared_for_plots.tsv"
# SCZ (PGC3) file:
SCZ_FILE <- "/Users/guillermocomesanacimadevila/Desktop/PhD/GenomicSEM/Part1/Data/SZ/PGC3_SCZ_wave3.harmonised_to_AD.prepared_for_plots.tsv"

# --- helpers ---------------------------------------------------------------
with_alpha <- function(cols, a = 0.80) {
  vapply(cols, function(c) grDevices::adjustcolor(c, alpha.f = a), "", USE.NAMES = FALSE)
}

prep_for_qqman <- function(dt) {
  # Ensure required columns and types: CHR (int), BP (int), P (numeric), SNP (char)
  need <- c("CHR","BP","P")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))
  # Coerce
  dt[, CHR := as.integer(CHR)]
  dt[, BP  := as.integer(BP)]
  dt[, P   := as.numeric(P)]
  # Filter invalid
  dt <- dt[is.finite(P) & P > 0 & !is.na(CHR) & !is.na(BP)]
  # SNP fallback
  if (!"SNP" %in% names(dt)) {
    if ("variant_id" %in% names(dt)) {
      dt[, SNP := as.character(variant_id)]
    } else {
      dt[, SNP := paste0(CHR, ":", BP)]
    }
  } else {
    dt[, SNP := as.character(SNP)]
    dt[is.na(SNP) | SNP == "", SNP := paste0(CHR, ":", BP)]
  }
  # Keep standard columns
  as.data.frame(dt[, .(SNP, CHR, BP, P)])
}

# --- read data -------------------------------------------------------------
ad_raw  <- data.table::fread(AD_FILE, sep = "\t", header = TRUE, showProgress = FALSE)
scz_raw <- data.table::fread(SCZ_FILE, sep = "\t", header = TRUE, showProgress = FALSE)

d_top    <- prep_for_qqman(ad_raw)   # AD top
d_bottom <- prep_for_qqman(scz_raw)  # SCZ bottom

# --- plotting --------------------------------------------------------------
plot_two_manhattans <- function(
    d_top, d_bottom,
    out    = "~/Downloads/manhattan_AD_top_SCZ_bottom.real.pdf",
    width  = 7.5,
    height = 5.8
) {
  # palette + thresholds
  blues   <- c("#0B3C5D", "#1F5A89", "#2F7FBF", "#66ACDC")
  blues_a <- with_alpha(blues, 0.80)
  gw_p  <- 5e-8;  sug_p <- 1e-5
  gw  <- -log10(gw_p);  sug <- -log10(sug_p)
  
  # fixed y-limit
  ymax <- 35
  
  # device
  try(graphics.off(), silent = TRUE)
  out_path <- path.expand(out)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  pdf(out_path, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  
  # base style
  op <- par(
    mfrow = c(2, 1),
    mgp   = c(2.0, 0.6, 0),
    oma   = c(0, 0, 0, 0),
    xaxs  = "i",
    tcl   = -0.25,
    lwd   = 1.0,
    lend  = "butt",
    bty   = "l"
  )
  on.exit(par(op), add = TRUE)
  
  add_overlays <- function() {
    usr <- par("usr")
    ybreaks <- pretty(c(0, usr[4]), n = 7)
    abline(h = ybreaks, col = "#E0E0E0", lwd = 0.6)
    abline(h = gw,  col = "#7D3C98", lty = 2, lwd = 1.6)
    abline(h = sug, col = "#BFC5CA", lty = 2, lwd = 1.2)
  }
  
  # ---------- TOP panel (A): no chr labels ----------
  par(mar = c(0.7, 3.8, 0.9, 0.8))  # tight bottom to reduce inter-panel gap
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
    chrlabs        = rep("", length(unique(d_top$CHR))),
    main           = ""
  )
  add_overlays()
  usr <- par("usr")
  text(x = usr[1] - 0.3 * diff(usr[1:2]), y = ymax * 0.9,
       labels = "A", font = 2, cex = 1.2, xpd = NA)
  
  # ---------- BOTTOM panel (B): keep chr labels ----------
  par(mar = c(2.6, 3.8, 0.5, 0.8))  # tight top to reduce gap
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

outfile <- plot_two_manhattans(
  d_top, d_bottom,
  out = "~/Downloads/manhattan_AD_top_SCZ_bottom.real.pdf"
)

cat("✔ Saved:", outfile, "\n")
