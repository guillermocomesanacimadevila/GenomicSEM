#!/usr/bin/env python3
# three_panel_manhattan_aesthetic.py
# 3 vertically stacked Manhattan plots with:
# - tighter vertical spacing
# - per-chromosome RAINBOW colours in a muted, publication-style palette
# - genome-wide significant points as triangles (same chromosome colour)
# - high default DPI (600)

import re, math, argparse, sys, gzip, io, csv
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import hsv_to_rgb, to_hex

# --- style ---
mpl.rcParams.update({
    "figure.dpi": 220,           # on-screen
    "savefig.dpi": 600,          # high-res export
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 13,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "axes.titlepad": 8,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
})
mpl.rcParams["pdf.compression"] = 0

# --- colours ---
# Build a muted "publication-tone" rainbow for chromosomes 1..22, X, Y, MT (25 total).
# Lower saturation and value than the previous version; slight lightness ramp to aid distinction.
def _rainbow_chr_palette(n_needed=25, base_s=0.55, base_v=0.82, v_ramp=0.08):
    """
    n_needed: number of colours required (25 for 1..22,X,Y,MT)
    base_s:   base saturation (0..1), lower = more muted
    base_v:   base value/brightness (0..1), lower = more muted
    v_ramp:   adds a small brightness ramp across hues for subtle contrast
    """
    hues = np.linspace(0.0, 1.0, n_needed, endpoint=False)
    # Apply a small sinusoidal lightness modulation to avoid adjacent hues looking identical in print
    phase = np.linspace(0, 2*np.pi, n_needed, endpoint=False)
    v = np.clip(base_v + v_ramp * 0.5 * np.sin(phase), 0.0, 1.0)
    s = np.full(n_needed, base_s)
    hsv = np.stack([hues, s, v], axis=1)
    rgb = hsv_to_rgb(hsv)
    return np.array([to_hex(c) for c in rgb])

CHR_COLOURS = _rainbow_chr_palette()  # 0->chr1, ..., 21->chr22, 22->X, 23->Y, 24->MT
COL_RED = "#B23A3A"   # slightly muted red for GW line
COL_BLUE = "#5E79C6"  # slightly muted blue for suggestive line

# --- io helpers ---
def _open_auto(path: str):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return io.StringIO(fh.read())
    return path

def _bump_csv_field_limit():
    try:
        csv.field_size_limit(sys.maxsize)
    except OverflowError:
        csv.field_size_limit(int(1e9))

# --- load + prep ---
def load_gwas_table(path: str) -> pd.DataFrame:
    _bump_csv_field_limit()
    src = _open_auto(path)
    try:
        df = pd.read_table(src, sep=None, engine="python", comment="#", dtype=str, low_memory=False)
    except Exception:
        try:
            df = pd.read_csv(src, sep="\t", engine="c", comment="#", dtype=str, low_memory=False)
        except Exception:
            df = pd.read_csv(src, sep=",", engine="c", comment="#", dtype=str, low_memory=False)

    colmap = {c.lower(): c for c in df.columns}
    def pick(*cands):
        for c in cands:
            if c in colmap: return colmap[c]
        return None

    chr_col = pick("chr", "chrom", "chromosome", "chrom_id")
    bp_col  = pick("bp", "pos", "position", "basepair", "base_pair")
    p_col   = pick("p", "pval", "pvalue", "p.value", "p_value", "p-val", "p-value")

    missing = [nm for nm,val in zip(["CHR","BP","P"], [chr_col,bp_col,p_col]) if val is None]
    if missing:
        raise ValueError(
            f"Missing required columns {missing} in file: {path}. "
            f"Found columns (first 12): {list(df.columns)[:12]}{'...' if len(df.columns)>12 else ''}"
        )

    out = pd.DataFrame({
        "CHR": df[chr_col],
        "BP":  df[bp_col],
        "P":   df[p_col],
    })
    for opt, aliases in {
        "SNP": ("snp","rsid","marker","id","markername"),
        "Effect": ("effect_allele","ea","a1","allele1","effect"),
        "Non_Effect": ("non_effect_allele","nea","a2","allele2","non_effect"),
        "is_snp": ("is_snp","snp_flag"),
        "is_ambiguous": ("is_ambiguous","ambiguous"),
    }.items():
        col = next((colmap[a] for a in aliases if a in colmap), None)
        if col is not None:
            out[opt] = df[col]
    return out

def prepare_for_plots(df: pd.DataFrame):
    df = df.copy()

    def _norm_chr(x):
        if pd.isna(x): return np.nan
        s = str(x).strip()
        s = re.sub(r"^chr", "", s, flags=re.IGNORECASE).upper()
        if s == "X": return 23
        if s == "Y": return 24
        if s in {"MT","M","MITO"}: return 25
        try: return int(float(s))
        except: return np.nan

    df["CHR_int"] = df["CHR"].map(_norm_chr).astype("Int64")
    df["BP_int"]  = pd.to_numeric(df["BP"], errors="coerce").astype("Int64")

    p = pd.to_numeric(df["P"], errors="coerce")
    p = p.where((p > 0) & (p <= 1), np.nan).clip(lower=1e-300)
    df["minus_log10p"] = -np.log10(p)

    df = df.dropna(subset=["CHR_int","BP_int","minus_log10p"]).copy()
    df.sort_values(["CHR_int","BP_int"], inplace=True)

    autos = sorted([c for c in df["CHR_int"].unique() if pd.notna(c) and c <= 22])
    others = [c for c in [23,24,25] if c in df["CHR_int"].unique()]
    chrom_order = autos + others

    xticks, xlabels, running = [], [], 0
    pos_cum = np.empty(len(df), dtype=np.int64)
    for c in chrom_order:
        mask = df["CHR_int"] == c
        if not mask.any(): continue
        bp = df.loc[mask, "BP_int"]
        bp_min, bp_max = int(bp.min()), int(bp.max())
        pos_cum[mask] = bp + running
        # center each visible chromosome block:
        xticks.append(running + (bp_min + bp_max) / 2.0)
        xlabels.append(str(c if c <= 22 else {23:"X",24:"Y",25:"MT"}[c]))
        running += (bp_max - bp_min + 1)

    df["pos_cum"] = pos_cum
    return df, xticks, xlabels, {"chrom_order": chrom_order}

def _to_bool_mask(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    s = series.astype(str).str.strip().str.lower()
    mask = s.isin(["1","true","t","yes","y"])
    mask = mask & ~s.isin(["nan","none","", "na"])
    return mask.fillna(False)

# --- plotting ---
def _chr_to_color(chr_int_series: pd.Series) -> np.ndarray:
    # Map 1..22->0..21, X->22, Y->23, MT->24 (wrap if needed)
    arr = chr_int_series.to_numpy(dtype=int)
    idx = np.where(arr <= 22, arr - 1,
                   np.where(arr == 23, 22,
                            np.where(arr == 24, 23, 24)))
    return CHR_COLOURS[idx % len(CHR_COLOURS)]

def manhattan_panel(
    ax, df_prep, xticks, xlabels,
    *, genomewide=5e-8, suggestive=1e-5,
    snp_only=True, drop_ambiguous=False,
    point_size=4.2, alpha=0.9,
    forced_ymax=None
):
    d = df_prep.copy()

    if snp_only and "is_snp" in d.columns:
        d = d[_to_bool_mask(d["is_snp"])]
    if drop_ambiguous and "is_ambiguous" in d.columns:
        d = d[~_to_bool_mask(d["is_ambiguous"])]

    # per-chromosome colours
    colours = _chr_to_color(d["CHR_int"])

    # base points
    ax.scatter(d["pos_cum"], d["minus_log10p"],
               s=point_size, c=colours, alpha=alpha,
               linewidths=0, rasterized=True, zorder=1)

    gw_y = -math.log10(genomewide)
    sug_y = -math.log10(suggestive)

    # significant points as triangles in their chromosome colours
    sig_mask = d["minus_log10p"] >= gw_y
    if np.any(sig_mask):
        ax.scatter(d.loc[sig_mask, "pos_cum"],
                   d.loc[sig_mask, "minus_log10p"],
                   s=max(point_size*2.0, 12),
                   marker="^",
                   c=colours[sig_mask],
                   edgecolor="white",
                   linewidths=0.35,
                   alpha=0.95,
                   zorder=3)

    # threshold lines
    ax.axhline(sug_y, color=COL_BLUE, linewidth=0.9, zorder=2)
    ax.axhline(gw_y, color=COL_RED,  linewidth=1.0, zorder=2)

    ymax = forced_ymax if forced_ymax is not None else float(
        max(np.nanpercentile(d["minus_log10p"], 99.9), 12.0) * 1.05
    )
    ax.set_ylim(0, ymax)
    ax.set_xlim(float(d["pos_cum"].min()), float(d["pos_cum"].max()))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.margins(x=0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.tick_params(direction="out", length=3.2, width=0.8)
    ax.set_rasterization_zorder(2)

def add_panel_label(ax, label):
    ax.text(0.01, 0.99, label, transform=ax.transAxes,
            ha="left", va="top", fontsize=14, fontweight="bold")

# --- 3-panel figure ---
def figure_three_manhattans_vertical(
    df1_prep, xt1, xl1,
    df2_prep, xt2, xl2,
    df3_prep, xt3, xl3,
    *,
    genomewide=5e-8, suggestive=1e-5,
    figsize=(12, 11.5),
    dpi=600,
    out_path=None, out_formats=("png","pdf","svg"),
    snp_only=True, drop_ambiguous=False,
    ymax_top=None, ymax_mid=None, ymax_bottom=None
):
    # Tight vertical spacing via gridspec hspace
    fig, axes = plt.subplots(
        nrows=3, ncols=1, figsize=figsize, dpi=dpi,
        sharex=False, constrained_layout=False,
        gridspec_kw={"hspace": 0.06}  # smaller -> closer panels
    )

    # TOP (A)
    manhattan_panel(axes[0], df1_prep, xt1, xl1,
                    genomewide=genomewide, suggestive=suggestive,
                    point_size=4.2, alpha=0.9, forced_ymax=ymax_top,
                    snp_only=snp_only, drop_ambiguous=drop_ambiguous)
    add_panel_label(axes[0], "A")
    axes[0].set_xlabel("")
    axes[0].set_xticks([])
    axes[0].set_xticklabels([])

    # MID (B)
    manhattan_panel(axes[1], df2_prep, xt2, xl2,
                    genomewide=genomewide, suggestive=suggestive,
                    point_size=4.2, alpha=0.9, forced_ymax=ymax_mid,
                    snp_only=snp_only, drop_ambiguous=drop_ambiguous)
    add_panel_label(axes[1], "B")
    axes[1].set_xlabel("")
    axes[1].set_xticks([])
    axes[1].set_xticklabels([])

    # BOTTOM (C): show chromosome axis
    manhattan_panel(axes[2], df3_prep, xt3, xl3,
                    genomewide=genomewide, suggestive=suggestive,
                    point_size=4.2, alpha=0.9, forced_ymax=ymax_bottom,
                    snp_only=snp_only, drop_ambiguous=drop_ambiguous)
    add_panel_label(axes[2], "C")
    axes[2].set_xlabel("Chromosome")

    if out_path:
        meta = {"Title": "Three-trait Manhattan", "Creator": "matplotlib"}
        for ext in out_formats:
            plt.savefig(f"{out_path}.{ext}", metadata=meta, dpi=dpi)

    return fig

# --- cli ---
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Three-panel Manhattan plot from three GWAS files."
    )
    p.add_argument("--trait1", required=True, help="Path to GWAS table for panel A (can be .gz)")
    p.add_argument("--trait2", required=True, help="Path to GWAS table for panel B (can be .gz)")
    p.add_argument("--trait3", required=True, help="Path to GWAS table for panel C (can be .gz)")
    p.add_argument("--out", default=None, help="Output path prefix (no extension). If omitted, just shows the plot.")
    p.add_argument("--genomewide", type=float, default=5e-8, help="Genome-wide threshold (default 5e-8)")
    p.add_argument("--suggestive", type=float, default=1e-5, help="Suggestive threshold (default 1e-5)")
    p.add_argument("--width", type=float, default=12.0, help="Figure width in inches")
    p.add_argument("--height", type=float, default=11.5, help="Figure height in inches")
    p.add_argument("--dpi", type=float, default=600.0, help="Figure and export DPI (default 600)")
    p.add_argument("--ymax-top", type=float, default=None, help="Force fixed y-axis max for TOP panel (A)")
    p.add_argument("--ymax-mid", type=float, default=None, help="Force fixed y-axis max for MIDDLE panel (B)")
    p.add_argument("--ymax-bottom", type=float, default=None, help="Force fixed y-axis max for BOTTOM panel (C)")
    p.add_argument("--no-snp-only", action="store_true", help="Do NOT restrict to rows flagged as SNPs when 'is_snp' present")
    p.add_argument("--drop-ambiguous", action="store_true", help="Drop rows flagged ambiguous if 'is_ambiguous' present")
    p.add_argument("--formats", type=str, default="png,pdf,svg",
                   help="Comma-separated export formats (e.g., 'png,svg' to skip PDF).")
    return p.parse_args(argv)

# --- main ---
if __name__ == "__main__":
    args = parse_args()

    df1_raw = load_gwas_table(args.trait1)
    df2_raw = load_gwas_table(args.trait2)
    df3_raw = load_gwas_table(args.trait3)

    prep1, xt1, xl1, _ = prepare_for_plots(df1_raw)
    prep2, xt2, xl2, _ = prepare_for_plots(df2_raw)
    prep3, xt3, xl3, _ = prepare_for_plots(df3_raw)

    formats = tuple([f.strip().lower() for f in args.formats.split(",") if f.strip()])

    figure_three_manhattans_vertical(
        prep1, xt1, xl1,
        prep2, xt2, xl2,
        prep3, xt3, xl3,
        genomewide=args.genomewide, suggestive=args.suggestive,
        figsize=(args.width, args.height),
        out_path=args.out,
        out_formats=formats,
        snp_only=not args.no_snp_only,
        drop_ambiguous=args.drop_ambiguous,
        ymax_top=args.ymax_top,
        ymax_mid=args.ymax_mid,
        ymax_bottom=args.ymax_bottom,
        dpi=int(args.dpi),
    )

    if args.out is None:
        plt.show()
