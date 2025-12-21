#!/usr/bin/env python3
import argparse
import pandas as pd
import subprocess
from pathlib import Path

# grab output of conjFDR
# grab snp hits,
# two dataframes (one AD, SCZ, - triple overlap = + LON)
# so 2 dfs for each snp based on pairwise traits
# for snp in sumstats_trait1:
#   clump r2 >= 0.6 per 1000Kb
# define lead SNP of locus
# around lead snp keep for each trait
# r2 >= 0.1 with lead over +/- 250Kb

def read_tsv(p):
    return pd.read_csv(p, sep='\t')

def write_tsv(df, p):
    df.to_csv(p, sep="\t", index=False)

def ensure_cols(df, cols):
    for col in cols:
        if col not in df.columns:
            raise ValueError(c)

def make_score_pairwise(df, score_col):
    ensure_cols(df, ["SNP", score_col])
    out = df[["SNP", score_col]].copy()
    out.columns = ["SNP", "P"]
    return out

def make_score_triple(df, c1, c2, how):
    ensure_cols(df, ["SNP", c1, c2])
    x = df[[c1, c2]].copy()
    if how == "max":
        p = x.max(axis=1)
    else:
        p = x.min(axis=1)
    out = pd.DataFrame({"SNP": df["SNP"], "P": p})
    return out

def run(cmd, cwd=None):
    subprocess.run(cmd, check=True, cwd=cwd)

def plink_clump(plink, ref_bfile, clump_in, out_prefix, r2, kb):
    cmd = [
        plink,
        "--bfile", ref_bfile,
        "--clump", str(clump_in),
        "--clump-field", "P",
        "--clump-p1", "1",
        "--clump-p2", "1",
        "--clump-r2", str(r2),
        "--clump-kb", str(kb),
        "--out", str(out_prefix),
    ]
    run(cmd)

def read_leads(clumped_path):
    df = pd.read_csv(clumped_path, sep=r"\s+")
    ensure_cols(df, ["SNP", "CHR", "BP"])
    df = df[["SNP", "CHR", "BP"]].dropna().drop_duplicates()
    df["locus_id"] = range(len(df))
    return df

def plink_ld(plink, ref_bfile, lead_snp, out_prefix, kb):
    cmd = [
        plink,
        "--bfile", ref_bfile,
        "--r2",
        "--ld-snp", lead_snp,
        "--ld-window-kb", str(kb),
        "--ld-window", "999999",
        "--ld-window-r2", "0",
        "--out", str(out_prefix),
    ]
    run(cmd)

def parse_ld_file(ld_path, r2_min):
    df = pd.read_csv(ld_path, sep=r"\s+")
    ensure_cols(df, ["SNP_A", "SNP_B", "CHR_A", "BP_A", "CHR_B", "BP_B", "R2"])
    df = df[df["R2"] >= r2_min].copy()
    if df.empty:
        return None, None
    snps = pd.concat([df["SNP_A"], df["SNP_B"]]).dropna().drop_duplicates()
    chr_vals = pd.concat([df["CHR_A"], df["CHR_B"]]).dropna().unique()
    chr_val = chr_vals[0]
    bp_min = int(min(df["BP_A"].min(), df["BP_B"].min()))
    bp_max = int(max(df["BP_A"].max(), df["BP_B"].max()))
    return snps, (chr_val, bp_min, bp_max)


def extract_sumstats(sumstats_path, snps):
    df = read_tsv(sumstats_path)
    ensure_cols(df, ["SNP"])
    keep = set(snps.tolist())
    out = df[df["SNP"].isin(keep)].copy()
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pairwise", "triple"], required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--plink", default="plink")
    ap.add_argument("--ref-bfile", required=True)
    ap.add_argument("--clump-r2", type=float, default=0.6)
    ap.add_argument("--clump-kb", type=int, default=1000)
    ap.add_argument("--ld-r2", type=float, default=0.1)
    ap.add_argument("--ld-kb", type=int, default=250)
    ap.add_argument("--score-col", default="conj_fdr")
    ap.add_argument("--triple-col1", default="conj_fdr_AD_LON")
    ap.add_argument("--triple-col2", default="conj_fdr_SCZ_LON")
    ap.add_argument("--triple-how", choices=["max", "min"], default="max")
    ap.add_argument("--sumstats", nargs="*", default=[])
    ap.add_argument("--sumstats-names", nargs="*", default=[])
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = read_tsv(args.hits).drop_duplicates(subset=["SNP"]).copy()
    if args.mode == "pairwise":
        clump_df = make_score_pairwise(hits, args.score_col)
    else:
        clump_df = make_score_triple(hits, args.triple_col1, args.triple_col2, args.triple_how)

    clump_in = out_dir / "clump_input.tsv"
    write_tsv(clump_df, clump_in)
    clump_out = out_dir / "clump"
    plink_clump(args.plink, args.ref_bfile, clump_in, clump_out, args.clump_r2, args.clump_kb)
    clumped_path = Path(str(clump_out) + ".clumped")
    leads = read_leads(clumped_path)
    write_tsv(leads, out_dir / "lead_snps.tsv")
    loci_rows = []
    for _, r in leads.iterrows():
        locus_id = int(r["locus_id"])
        lead = r["SNP"]
        pref = out_dir / f"locus_{locus_id}"
        plink_ld(args.plink, args.ref_bfile, lead, pref, args.ld_kb)
        ld_path = Path(str(pref) + ".ld")
        snps, coords = parse_ld_file(ld_path, args.ld_r2)
        if snps is None:
            continue
        chr_val, start, end = coords
        loci_rows.append({"locus_id": locus_id, "lead_snp": lead, "chr": chr_val, "start": start, "end": end, "n_snps": len(snps)})
        write_tsv(pd.DataFrame({"SNP": snps}), out_dir / f"locus_{locus_id}_snps.tsv")
        if len(args.sumstats) == len(args.sumstats_names) and len(args.sumstats) > 0:
            for p, name in zip(args.sumstats, args.sumstats_names):
                loc = extract_sumstats(p, snps)
                write_tsv(loc, out_dir / f"locus_{locus_id}_{name}.tsv")

    if len(loci_rows) > 0:
        write_tsv(pd.DataFrame(loci_rows), out_dir / "locus_coords.tsv")
    else:
        write_tsv(pd.DataFrame(columns=["locus_id","lead_snp","chr","start","end","n_snps"]), out_dir / "locus_coords.tsv")

if __name__ == "__main__":
    main()