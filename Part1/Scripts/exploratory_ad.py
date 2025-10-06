#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def convert(path: str) -> str:
    df = pd.read_csv(path, sep=r"\s+")
    tsv_path = os.path.splitext(path)[0] + ".tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path

def qc_checks(df: pd.DataFrame, drop_indels: bool = False) -> pd.DataFrame:
    print("===== QC Summary =====")
    print(f"Number of SNPs: {df.shape[0]:,}")
    print(f"Number of columns: {df.shape[1]}")
    print("\nColumn names:", list(df.columns))
    print("\nMissing values per column:")
    print(df.isna().sum())

    cols_lower = [c.lower() for c in df.columns]
    if "beta" in cols_lower:
        print("\nEffect size column: BETA")
    elif "or" in cols_lower:
        print("\nEffect size column: OR")
    else:
        print("\nNo BETA/OR column found")

    pcol = "P" if "P" in df.columns else None
    if pcol:
        vmax = pd.to_numeric(df[pcol], errors="coerce").max()
        print(f"{pcol} column looks like {'-log10(P)' if (pd.notna(vmax) and vmax > 1) else 'raw P-values'}")
    else:
        print("No P/PVAL column found")

    allele_pairs = [
        ("Effect", "Non_Effect"),
        ("A1", "A2"),
        ("EA", "NEA"),
        ("ALT", "REF"),
        ("ALLELE1", "ALLELE2"),
    ]
    a1, a2 = None, None
    for c1, c2 in allele_pairs:
        if c1 in df.columns and c2 in df.columns:
            a1, a2 = c1, c2
            break

    if a1 and a2:
        bases = {"A", "C", "G", "T"}
        a1s = df[a1].astype(str).str.upper()
        a2s = df[a2].astype(str).str.upper()
        is_single = a1s.isin(bases) & a2s.isin(bases) & (a1s.str.len() == 1) & (a2s.str.len() == 1)
        has_gap = a1s.str.contains("-", na=False) | a2s.str.contains("-", na=False)
        is_snp = is_single & (~has_gap)
        is_indel = ~is_snp
        n_indels = int(is_indel.sum())
        n_snps_before = int(is_snp.sum())
        print(f"\nINDELs detected: {n_indels:,}")
        print(f"SNPs (before drop): {n_snps_before:,}")
        if drop_indels:
            df = df[is_snp].copy()
            print(f"SNPs (after drop): {df.shape[0]:,}")
    else:
        print("\nINDEL check skipped: no allele columns found")

    print("=======================")
    return df

def prepare_for_plots(df: pd.DataFrame):
    df_out = df.copy()

    def _norm_chr(x):
        if pd.isna(x): return np.nan
        s = str(x).strip()
        s = re.sub(r"^chr", "", s, flags=re.IGNORECASE).upper()
        if s == "X": return 23
        if s == "Y": return 24
        if s in {"MT", "M", "MITO"}: return 25
        try:
            return int(float(s))
        except Exception:
            return np.nan

    df_out["CHR_int"] = df_out["CHR"].map(_norm_chr).astype("Int64")
    df_out["BP_int"] = pd.to_numeric(df_out["BP"], errors="coerce").astype("Int64")

    pvals = pd.to_numeric(df_out["P"], errors="coerce")
    tiny = 1e-300
    pvals = pvals.where((pvals > 0) & (pvals <= 1), np.nan).clip(lower=tiny)
    df_out["minus_log10p"] = -np.log10(pvals)

    eff = df_out["Effect"].astype(str)
    non = df_out["Non_Effect"].astype(str)
    bases = {"A", "C", "G", "T"}

    is_snp = eff.isin(bases) & non.isin(bases) & (eff.str.len() == 1) & (non.str.len() == 1)
    is_indel = ~is_snp | eff.str.contains("-", na=False) | non.str.contains("-", na=False)

    ambiguous = (is_snp & (
        ((eff == "A") & (non == "T")) |
        ((eff == "T") & (non == "A")) |
        ((eff == "C") & (non == "G")) |
        ((eff == "G") & (non == "C"))
    ))

    df_out["is_snp"] = is_snp.fillna(False)
    df_out["is_indel"] = is_indel.fillna(False)
    df_out["is_ambiguous"] = ambiguous.fillna(False)

    snp_id = df_out["SNP"].astype(str)
    rs_mask = snp_id.str.startswith("rs", na=False)
    fallback = (
        df_out["CHR_int"].astype(str)
        + ":" + df_out["BP_int"].astype(str)
        + ":" + non + "_" + eff
    )
    df_out["variant_id"] = np.where(
        rs_mask, snp_id,
        np.where(df_out["Test_MarkerName"].notna(),
                 df_out["Test_MarkerName"].astype(str),
                 fallback)
    )

    if "Effect_allele_freq" in df_out.columns:
        eaf = pd.to_numeric(df_out["Effect_allele_freq"], errors="coerce")
        df_out["eaf"] = eaf
        df_out["maf"] = np.minimum(eaf, 1 - eaf)
    else:
        df_out["eaf"] = np.nan
        df_out["maf"] = np.nan

    before = len(df_out)
    df_out = df_out.dropna(subset=["CHR_int", "BP_int", "minus_log10p"]).copy()
    n_dropped = before - len(df_out)

    autos = sorted([c for c in df_out["CHR_int"].unique() if pd.notna(c) and c <= 22])
    others = [c for c in [23, 24, 25] if c in df_out["CHR_int"].unique()]
    chrom_order = autos + others

    df_out.sort_values(["CHR_int", "BP_int"], inplace=True)
    xticks, xlabels = [], []
    running = 0
    for c in chrom_order:
        mask = df_out["CHR_int"] == c
        if not mask.any(): continue
        bp_min, bp_max = df_out.loc[mask, "BP_int"].agg(["min", "max"])
        df_out.loc[mask, "pos_cum"] = df_out.loc[mask, "BP_int"] + running
        mid = running + (bp_max - bp_min) / 2
        xticks.append(mid)
        xlabels.append(str(c if c <= 22 else {23: "X", 24: "Y", 25: "MT"}[c]))
        running += (bp_max - bp_min + 1)

    df_out["pos_cum"] = df_out["pos_cum"].astype(np.int64)

    meta = {"n_dropped": n_dropped, "chrom_order": chrom_order}
    return df_out, xticks, xlabels, meta

def manhattan_plot(
    df_prep: pd.DataFrame,
    xticks: list,
    xlabels: list,
    *,
    genomewide: float = 5e-8,
    suggestive: float = 1e-5,
    snp_only: bool = False,
    drop_ambiguous: bool = False,
    title: str = "Manhattan plot",
    point_size: float = 1.0,
    alpha: float = 0.9,
    ymax: float | None = None,
    figsize: tuple = (12, 5),
    dpi: int = 300,
    out_path: str | None = None,          
    out_formats: tuple = ("png", "pdf"),  
    annotate_top_n: int | None = None,    
    annotate_threshold: float | None = None,  
    highlight: list[str] | None = None,  
    font_size_axis: int = 10,
    font_size_tick: int = 9,
    font_size_title: int = 12
):

    d = df_prep.copy()
    if snp_only and "is_snp" in d:
        d = d[d["is_snp"]]
    if drop_ambiguous and "is_ambiguous" in d:
        d = d[~d["is_ambiguous"]]

    y = d["minus_log10p"].to_numpy()
    if ymax is None:
        ymax = max(np.nanpercentile(y, 99.9), 40.0) * 1.05

    gw_line = -np.log10(genomewide) if genomewide else None
    sug_line = -np.log10(suggestive) if suggestive else None

    base_colors = np.array([
        "#4C78A8", "#9EB7D4"  
    ])
    chr_vals = d["CHR_int"].to_numpy()
    color_idx = (chr_vals % 2).astype(int)
    colors = base_colors[color_idx]

    highlight_set = set(highlight or [])
    if "variant_id" in d.columns and highlight_set:
        is_high = d["variant_id"].isin(highlight_set)
    else:
        is_high = pd.Series(False, index=d.index)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    if len(xlabels) == len(xticks) and len(xticks) > 1:
        mids = np.array(xticks, dtype=float)
        edges = np.empty(len(mids) + 1, dtype=float)
        edges[1:-1] = (mids[1:] + mids[:-1]) / 2
        edges[0] = float(d["pos_cum"].min())
        edges[-1] = float(d["pos_cum"].max())
        for i in range(len(mids)):
            if i % 2 == 0:
                ax.axvspan(edges[i], edges[i+1], facecolor="#F5F7FA", alpha=0.9, zorder=0)

    ax.scatter(
        d["pos_cum"],
        d["minus_log10p"],
        s=point_size,
        alpha=alpha,
        c=colors,
        linewidths=0,
        rasterized=True
    )

    if is_high.any():
        ax.scatter(
            d.loc[is_high, "pos_cum"],
            d.loc[is_high, "minus_log10p"],
            s=max(point_size * 2, 6),
            alpha=1.0,
            c="#D62728",    
            edgecolor="white",
            linewidths=0.25,
            rasterized=True,
            zorder=3
        )

    def _hline(yval, ls, col, label):
        ax.axhline(yval, linestyle=ls, linewidth=1.0, color=col, alpha=0.9)
        ax.text(
            x=ax.get_xlim()[1], y=yval,
            s=label,
            va="bottom", ha="right",
            fontsize=font_size_tick, color=col,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6)
        )

    if sug_line:
        _hline(sug_line, "--", "#7F7F7F", f"Suggestive (p={suggestive:g})")
    if gw_line:
        _hline(gw_line, "-", "#C00000", f"Genome-wide (p={genomewide:g})")

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=font_size_tick)
    ax.set_xlim(float(d["pos_cum"].min()), float(d["pos_cum"].max()))
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Chromosome", fontsize=font_size_axis)
    ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=font_size_axis)
    ax.set_title(title, fontsize=font_size_title, pad=8)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.margins(x=0)

    if annotate_top_n or annotate_threshold:
        cand = d.copy()
        if annotate_threshold is not None:
            cand = cand[cand["minus_log10p"] >= float(annotate_threshold)]
        cand = cand.sort_values("minus_log10p", ascending=False)
        if annotate_top_n is not None:
            cand = cand.head(int(annotate_top_n))

        if "variant_id" in cand.columns:
            for _, row in cand.iterrows():
                ax.annotate(
                    text=str(row["variant_id"]),
                    xy=(float(row["pos_cum"]), float(row["minus_log10p"])),
                    xytext=(0, 6),
                    textcoords="offset points",
                    fontsize=max(font_size_tick - 1, 8),
                    ha="center",
                    va="bottom",
                    color="#111111",
                    clip_on=True
                )

    fig.tight_layout()

    if out_path:
        for ext in out_formats:
            fig.savefig(
                f"{out_path}.{ext}",
                dpi=dpi,
                bbox_inches="tight",
                pad_inches=0.05,
                transparent=False
            )

    return fig, ax


if __name__ == "__main__":
    path_txt = "/Users/guillermocomesanacimadevila/Desktop/PhD/Part1/Data/AD/Kunkle_etal_2019_IGAP_Summary_statistics_published.txt"
    df = pd.read_csv(path_txt, sep=r"\s+")
    df = qc_checks(df, drop_indels=True)

    print("\n[1/3] Writing RAW TSV (no -log10 transform)...")
    raw_tsv_path = convert(path_txt)
    print(f"RAW TSV: {raw_tsv_path}")

    print("\n[2/3] Preparing DataFrame for plots...")
    df_prep, xticks, xlabels, meta = prepare_for_plots(df)
    print("Prepared rows:", df_prep.shape[0], "| Dropped:", meta["n_dropped"])
    print("Chromosomes detected:", meta["chrom_order"])

    base_no_ext, _ = os.path.splitext(path_txt)
    prepared_tsv_path = base_no_ext + ".prepared_for_plots.tsv"
    df_prep.to_csv(prepared_tsv_path, sep="\t", index=False)
    print(f"PREPARED TSV: {prepared_tsv_path}")

    print("\n[3/3] Rendering Manhattan plot...")
    #manhattan_plot(
        #df_prep, xticks, xlabels,
        #genomewide=5e-8, suggestive=1e-5,
        #snp_only=True,
        #drop_ambiguous=False,
        #title="IGAP 2019 (SNP-only)",
        #point_size=0.6,
        #out_path="manhattan_IGAP2019",
        #out_formats=("png", "pdf"),
        #annotate_top_n=10,
        #annotate_threshold=7.3
    #)
    print("Done.")



# TO DO´s
# Do all of this but for the SZ GWAS
# Use Genomic SEM package for LDSC (R)
# Get genetic correlation between both GWAS - get a small table + heatmap 22
# r(theta) should be 0 or near 0 
# Also make sure Beta is not -log10(p) when slapping data into GenomicSEM LDSC