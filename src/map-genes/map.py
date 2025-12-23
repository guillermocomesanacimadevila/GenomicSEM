#!/usr/bin/env python3
import pandas as pd
import numpy as np
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

def map_ensembl_to_symbol(df1: pd.DataFrame, ref: pd.DataFrame, df1_col: str, ensembl_col: str, symbol_col:str) -> pd.DataFrame:
    df1 = df1.copy()
    ref = ref.copy()
    ref = reformat_ensembl_file(ref)
    df1["symbol"] = df1[df1_col].map(
        ref.set_index(ensembl_col)[symbol_col]
    )
    return df1

# =========
# =========

def extract_genes_pairwise(df1: pd.DataFrame, df2: pd.DataFrame, gene_col:str, pheno1_id:str, pheno2_id:str):
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
    print(f"Raw number of shared genes in this locus: {df[gene_col].unique()}")
    return df

def extract_genes_triple_overlap(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    df3: pd.DataFrame,
    gene_col,
    pheno1_id: str,
    pheno2_id: str,
    pheno3_id: str,
):
    df1 = df1.copy()
    df2 = df2.copy()
    df3 = df3.copy()

    def clean_gene(df: pd.DataFrame) -> pd.DataFrame:
        m = df[gene_col].astype(str).str.strip().str.lower().ne("nan")
        return df[m]

    df1 = clean_gene(df1)
    df2 = clean_gene(df2)
    df3 = clean_gene(df3)
    df12 = df1.merge(
        df2,
        on=gene_col,
        how="inner",
        suffixes=(f"_{pheno1_id}", f"_{pheno2_id}")
    )

    df3_renamed = df3.rename(
        columns={c: f"{c}_{pheno3_id}" for c in df3.columns if c != gene_col}
    )

    df123 = df12.merge(
        df3_renamed,
        on=gene_col,
        how="inner"
    )

    df123[gene_col] = df123[gene_col].replace({"FAM63B": "MINDY1", "FAM63A": "MINDY2"})
    print(f"Raw number of shared genes in this locus: {df123[gene_col].unique()}")
    return df123



# return dataframe and then call funct in main() -> .to_csv(sep="\t")
# ensure cols == harmonised to each pheno
def map_pairwise_genes_to_magma():
    return

# ensure cols == harmonised to each pheno (1-3)
def map_triple_overlap_genes_to_magma():
    return


# make if statement to check whether how many traits
# if 2 = do it for pairwise
# if 3 = triple overlap
def ensure_significance():
    return


# FRD - != avoid LD confounding
def multiple_testing():
    return








PATH = Path("/Users/c24102394/Desktop/PhD/DiscoveryPipeline/outputs/fuma_post")
ENSEMBL = parse_df("/Users/c24102394/ensemble/mart_export.txt")

ad0 = parse_df(PATH / "AD/locus_0/snps.txt")
scz0 = parse_df(PATH / "SCZ/locus_0/snps.txt")

ad_eqtl = parse_df(PATH / "AD/locus_0/eqtl.txt")
scz_eqtl = parse_df(PATH / "SCZ/locus_0/eqtl.txt")

ad_ci = parse_df(PATH / "AD/locus_0/ci.txt")
scz_ci = parse_df(PATH / "SCZ/locus_0/ci.txt")


df = extract_genes_pairwise(ad0, scz0, "nearestGene", "AD", "SCZ")
df2 = extract_genes_pairwise(ad_eqtl, scz_eqtl, "symbol", "AD", "SCZ")

df3_ad = map_ensembl_to_symbol(ad_ci, ENSEMBL, "genes", "ensembl_gene_id", "gene_symbol")
df3_scz = map_ensembl_to_symbol(scz_ci, ENSEMBL, "genes", "ensembl_gene_id", "gene_symbol")
df3 = extract_genes_pairwise(df3_ad, df3_scz, "symbol", "AD", "SCZ")


# triple overlap
ad2 = parse_df(PATH  / "AD/locus_2/snps.txt")
scz2 = parse_df(PATH / "SCZ/locus_2/snps.txt" )
lon2 = parse_df(PATH / "LON/locus_0/snps.txt")

ad2_eqtl = parse_df(PATH  / "AD/locus_2/eqtl.txt")
scz2_eqtl = parse_df(PATH / "SCZ/locus_2/eqtl.txt")
lon2_eqtl = parse_df(PATH / "LON/locus_0/eqtl.txt")

ad2_ci = parse_df(PATH  / "AD/locus_2/ci.txt")
scz2_ci = parse_df(PATH / "SCZ/locus_2/ci.txt")
lon2_ci = parse_df(PATH / "LON/locus_0/ci.txt")

triple_ci_ad = map_ensembl_to_symbol(ad2_ci, ENSEMBL, "genes", "ensembl_gene_id", "gene_symbol")
triple_ci_scz = map_ensembl_to_symbol(scz2_ci, ENSEMBL, "genes", "ensembl_gene_id", "gene_symbol")
triple_ci_lon = map_ensembl_to_symbol(lon2_ci, ENSEMBL, "genes", "ensembl_gene_id", "gene_symbol")

triple = extract_genes_triple_overlap(ad2, scz2, lon2, "nearestGene", "AD", "SCZ", "LON")
triple_eqtl = extract_genes_triple_overlap(ad2_eqtl, scz2_eqtl, lon2_eqtl, "symbol", "AD", "SCZ", "LON")
triple_ci = extract_genes_triple_overlap(triple_ci_ad, triple_ci_scz, triple_ci_lon, "symbol", "AD", "SCZ", "LON")

