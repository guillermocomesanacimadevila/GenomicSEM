#!/usr/bin/env python3

import re, math, argparse, sys, gzip, io, csv
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

mpl.rcParams.update({
    "figure.dpi": 220,
    "savefig.dpi": 300,
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

COL_BLACK = "#222222"
COL_GREY  = "#B7B7B7"
COL_SIG   = "#14B5FF"
COL_RED   = "#C00000"
COL_BLUE  = "#6A8DFF"

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
        "SNP": ("snp","rsid","marker","id"),
        "Effect": ("effect_allele","ea","a1","allele1"),
        "Non_Effect": ("non_effect_allele","nea","a2","allele2"),
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
        bp_min, bp_max = bp.min(), bp.max()
        pos_cum[mask] = bp + running
        xticks.append(running + (bp_max - bp_min)/2)
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

    chr_vals = d["CHR_int"].to_numpy()
    colors = np.where((chr_vals % 2) == 0, COL_BLACK, COL_GREY)

    ax.scatter(d["pos_cum"], d["minus_log10p"],
               s=point_size, c=colors, alpha=alpha,
               linewidths=0, rasterized=True, zorder=1)

    gw_y = -math.log10(genomewide)
    sug_y = -math.log10(suggestive)

    sig_mask = d["minus_log10p"] >= gw_y
    if np.any(sig_mask):
        ax.scatter(d.loc[sig_mask, "pos_cum"],
                   d.loc[sig_mask, "minus_log10p"],
                   s=max(point_size*2.0, 12),
                   marker="D", c=COL_SIG, edgecolor="white",
                   linewidths=0.35, alpha=0.95, zorder=3)

    ax.axhline(sug_y, color=COL_BLUE, linewidth=0.9, zorder=2)
    ax.axhline(gw_y, color=COL_RED,  linewidth=1.0, zorder=2)

    ymax = forced_ymax if forced_ymax is not None else float(max(np.nanpercentile(d["minus_log10p"], 99.9), 12.0) * 1.05)
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

def figure_two_manhattans_vertical(
    df1_prep, xt1, xl1,
    df2_prep, xt2, xl2,
    *,
    genomewide=5e-8, suggestive=1e-5,
    figsize=(12, 8),
    dpi=300,
    out_path=None, out_formats=("png","pdf","svg"),
    snp_only=True, drop_ambiguous=False,
    forced_ymax=None
):
    if forced_ymax is not None:
        ymax = float(forced_ymax)
    else:
        ymax = float(max(
            np.nanpercentile(np.r_[df1_prep["minus_log10p"], df2_prep["minus_log10p"]], 99.9),
            12.0
        ) * 1.05)

    fig, (ax1, ax2) = plt.subplots(
        nrows=2, ncols=1, figsize=figsize, dpi=dpi,
        sharex=False, constrained_layout=True
    )

    manhattan_panel(ax1, df1_prep, xt1, xl1,
                    genomewide=genomewide, suggestive=suggestive,
                    point_size=4.2, alpha=0.9, forced_ymax=ymax,
                    snp_only=snp_only, drop_ambiguous=drop_ambiguous)
    add_panel_label(ax1, "A")
    ax1.set_xlabel("")

    manhattan_panel(ax2, df2_prep, xt2, xl2,
                    genomewide=genomewide, suggestive=suggestive,
                    point_size=4.2, alpha=0.9, forced_ymax=ymax,
                    snp_only=snp_only, drop_ambiguous=drop_ambiguous)
    add_panel_label(ax2, "B")
    ax2.set_xlabel("Chromosome")

    if out_path:
        meta = {"Title": "Two-trait Manhattan", "Creator": "matplotlib"}
        for ext in out_formats:
            fig.savefig(f"{out_path}.{ext}", metadata=meta, dpi=dpi)

    return fig

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Two-panel Manhattan plot from two GWAS files."
    )
    p.add_argument("--trait1", required=True, help="Path to GWAS table for panel A (can be .gz)")
    p.add_argument("--trait2", required=True, help="Path to GWAS table for panel B (can be .gz)")
    p.add_argument("--out", default=None, help="Output path prefix (without extension). If omitted, just shows the plot.")
    p.add_argument("--genomewide", type=float, default=5e-8, help="Genome-wide threshold (default 5e-8)")
    p.add_argument("--suggestive", type=float, default=1e-5, help="Suggestive threshold (default 1e-5)")
    p.add_argument("--width", type=float, default=12.0, help="Figure width in inches")
    p.add_argument("--height", type=float, default=8.0, help="Figure height in inches")
    p.add_argument("--dpi", type=float, default=300.0, help="Figure and export DPI (default 300)")
    p.add_argument("--ymax", type=float, default=None, help="Force fixed y-axis max for both panels (optional)")
    p.add_argument("--no-snp-only", action="store_true", help="Do NOT restrict to rows flagged as SNPs when 'is_snp' present")
    p.add_argument("--drop-ambiguous", action="store_true", help="Drop rows flagged ambiguous if 'is_ambiguous' present")
    p.add_argument("--formats", type=str, default="png,pdf,svg",
                   help="Comma-separated export formats (e.g., 'png,svg' to skip PDF).")
    return p.parse_args(argv)

if __name__ == "__main__":
    args = parse_args()

    df1_raw = load_gwas_table(args.trait1)
    df2_raw = load_gwas_table(args.trait2)

    prep1, xt1, xl1, _ = prepare_for_plots(df1_raw)
    prep2, xt2, xl2, _ = prepare_for_plots(df2_raw)

    formats = tuple([f.strip().lower() for f in args.formats.split(",") if f.strip()])

    figure_two_manhattans_vertical(
        prep1, xt1, xl1,
        prep2, xt2, xl2,
        genomewide=args.genomewide, suggestive=args.suggestive,
        figsize=(args.width, args.height),
        out_path=args.out,
        out_formats=formats,
        snp_only=not args.no_snp_only,
        drop_ambiguous=args.drop_ambiguous,
        forced_ymax=args.ymax,
        dpi=args.dpi,
    )

    if args.out is None:
        plt.show()