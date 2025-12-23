#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse

# for nf
# -> mkdir ${workflow.launchDir}/outputs/gene-mappings
# for loc in {0..2}
# do
#   mkdir -p "${workflow.launchDir}/outputs/gene-mappings/${loc}"
# done
# files in dirs -> eQTL / HiC / positional -> mappings per locus

# args
# out_dir
# trait1_snps, trait1_eqtl, trait1_3dci
# trait2_snps, trait2_eqtl, trait2_3dci
# trait3_snps, trait3_eqtl, trait3_3dci
# magma_trait1, magma_trait2, magma_trait3
# e-magma_trait1, e-magma_trait2, e-magma_trait3
# h-magma_trait1, h-magma_trait2, h-magma_trait3

# Gene mapping per locus
# Criteria: (positional mapping - +/- 10Kb)
# * Gene set 1
# - Positional (within 10Kb) - Map all genes from SNP list of each locus
# - Positional - Ensure overlap between traits involved - if == overlap -> match
#
# * Gene set 2 (cis-eQTL mapping)
# - cis-eQTLs - FDR significant eQTL hits for region across all SNPs within locus
# - cis-eQTLs - Ensure eQTL overlaps across all traits involved - if == overlap -> match
#
# * Gene set 3 (3D-CI interaction mapping)
# - 3D-CI - FDR significant HiC hits for region across all SNPs within locus
# - 3D-CI - Ensure HiC overlaps across all traits involved - if == overlap -> match

def parse_df(path):
    return pd.read_csv(path, sep="\t")

# need to download the ensembl gene code with respective symbol
# BioMart - https://www.ensembl.org/biomart/martview/5710f6f7a825940b96468f0555d548e9
def reformat_ensembl_file(df):
    df = df.copy()
    df.columns = [
        "ensembl_gene_id",
        "ensembl_gene_id_version",
        "gene_symbol",
    ]
    return df

def map_ensembl_to_symbol(df1: pd.DataFrame, ref: pd.DataFrame, df1_col: str, ensembl_col: str, symbol_col: str) -> pd.DataFrame:
    df1 = df1.copy()
    ref = ref.copy()
    ref = reformat_ensembl_file(ref)
    df1["symbol"] = df1[df1_col].map(
        ref.set_index(ensembl_col)[symbol_col]
    )
    return df1

# =========
# =========

def extract_genes_pairwise(df1: pd.DataFrame, df2: pd.DataFrame, gene_col: str, pheno1_id: str, pheno2_id: str):
    df1 = df1.copy()
    df2 = df2.copy()
    df1 = df1[df1[gene_col].astype(str).str.strip().str.lower() != "nan"]
    df2 = df2[df2[gene_col].astype(str).str.strip().str.lower() != "nan"]
    df = df1.merge(
        df2,
        on=gene_col,
        how="inner",
        suffixes = (f"_{pheno1_id}", f"_{pheno2_id}")
    )
    df[gene_col] = df[gene_col].replace({"FAM63B": "MINDY1"})
    df[gene_col] = df[gene_col].replace({"FAM63A": "MINDY2"})
    print(f"[pairwise] shared genes ({gene_col}): {df[gene_col].nunique()}")
    return df

def extract_genes_triple_overlap(df1, df2, df3, gene_col, pheno1_id, pheno2_id, pheno3_id):
    df1 = df1.copy()
    df2 = df2.copy()
    df3 = df3.copy()

    def clean(df):
        s = df[gene_col].astype(str).str.strip()
        return df[s.ne("") & s.str.lower().ne("nan")]

    df1 = clean(df1)
    df2 = clean(df2)
    df3 = clean(df3)
    df12 = df1.merge(
        df2,
        on=gene_col,
        how="inner",
        suffixes=(f"_{pheno1_id}", f"_{pheno2_id}")
    )

    df3 = df3.rename(columns={c: f"{c}_{pheno3_id}" for c in df3.columns if c != gene_col})
    df = df12.merge(df3, on=gene_col, how="inner")
    df[gene_col] = df[gene_col].replace({"FAM63B": "MINDY1", "FAM63A": "MINDY2"})
    print(df[gene_col].nunique())
    return df

def save_outputs(df, gene_col, out_prefix):
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(f"{out_prefix}_overlap.tsv", sep="\t", index=False)
    genes = (
        df[gene_col]
        .astype(str)
        .str.strip()
        .replace({"nan": pd.NA, "": pd.NA})
        .dropna()
        .drop_duplicates()
        .sort_values()
    )
    genes.to_frame(gene_col).to_csv(f"{out_prefix}_genes.tsv", sep="\t", index=False)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno1_id", required=True)
    p.add_argument("--pheno2_id", required=True)
    p.add_argument("--pheno3_id")
    p.add_argument("--pheno1_snps", required=True)
    p.add_argument("--pheno1_eqtl", required=True)
    p.add_argument("--pheno1_ci", required=True)
    p.add_argument("--pheno2_snps", required=True)
    p.add_argument("--pheno2_eqtl", required=True)
    p.add_argument("--pheno2_ci", required=True)
    p.add_argument("--pheno3_snps")
    p.add_argument("--pheno3_eqtl")
    p.add_argument("--pheno3_ci")
    p.add_argument("--ensembl_ref", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--do_triple", action="store_true")
    args = p.parse_args()
    out = Path(args.out_dir)
    ensembl = parse_df(args.ensembl_ref)
    ad_snps = parse_df(args.pheno1_snps)
    scz_snps = parse_df(args.pheno2_snps)
    ad_eqtl = parse_df(args.pheno1_eqtl)
    scz_eqtl = parse_df(args.pheno2_eqtl)
    ad_ci = parse_df(args.pheno1_ci)
    scz_ci = parse_df(args.pheno2_ci)
    pos = extract_genes_pairwise(ad_snps, scz_snps, "nearestGene", args.pheno1_id, args.pheno2_id)
    save_outputs(pos, "nearestGene", out / f"{args.pheno1_id}__{args.pheno2_id}_positional")
    eqtl = extract_genes_pairwise(ad_eqtl, scz_eqtl, "symbol", args.pheno1_id, args.pheno2_id)
    save_outputs(eqtl, "symbol", out / f"{args.pheno1_id}__{args.pheno2_id}_eqtl")
    ad_ci = map_ensembl_to_symbol(ad_ci, ensembl, "genes", "ensembl_gene_id", "gene_symbol")
    scz_ci = map_ensembl_to_symbol(scz_ci, ensembl, "genes", "ensembl_gene_id", "gene_symbol")
    ci = extract_genes_pairwise(ad_ci, scz_ci, "symbol", args.pheno1_id, args.pheno2_id)
    save_outputs(ci, "symbol", out / f"{args.pheno1_id}__{args.pheno2_id}_ci")
    if args.do_triple:
        lon_snps = parse_df(args.pheno3_snps)
        lon_eqtl = parse_df(args.pheno3_eqtl)
        lon_ci = parse_df(args.pheno3_ci)
        pos3 = extract_genes_triple_overlap(
            ad_snps, scz_snps, lon_snps, "nearestGene",
            args.pheno1_id, args.pheno2_id, args.pheno3_id
        )
        save_outputs(pos3, "nearestGene", out / f"{args.pheno1_id}__{args.pheno2_id}__{args.pheno3_id}_positional")
        eqtl3 = extract_genes_triple_overlap(
            ad_eqtl, scz_eqtl, lon_eqtl, "symbol",
            args.pheno1_id, args.pheno2_id, args.pheno3_id
        )
        save_outputs(eqtl3, "symbol", out / f"{args.pheno1_id}__{args.pheno2_id}__{args.pheno3_id}_eqtl")
        lon_ci = map_ensembl_to_symbol(lon_ci, ensembl, "genes", "ensembl_gene_id", "gene_symbol")
        ci3 = extract_genes_triple_overlap(
            ad_ci, scz_ci, lon_ci, "symbol",
            args.pheno1_id, args.pheno2_id, args.pheno3_id
        )
        save_outputs(ci3, "symbol", out / f"{args.pheno1_id}__{args.pheno2_id}__{args.pheno3_id}_ci")

if __name__ == "__main__":
    main()
