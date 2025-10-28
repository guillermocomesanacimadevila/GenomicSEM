## Run MAGMA on current data

```bash
MAGMA
├── CHANGELOG
├── magma
├── magma_v1.10_mac.zip
├── manual_v1.10.pdf
├── README
└── ref
    ├── g1000_eur.bed
    ├── g1000_eur.bim
    ├── g1000_eur.fam
    ├── g1000_eur.synonyms
    ├── NCBI37.3.gene.loc
    ├── README
    ├── README2
    └── REPORT
```


### 1. Create MAGMA output dir

```bash
cd /Users/c24102394/Desktop/PhD/AD_SZ_genes/LAVA_LDSC_nature/outputs/lava
mkdir -p magma_runs && cd magma_runs
```


### 2. Create genes_b37.genes.annot by mapping 1kG EUR SNPs to NCBI37.3 gene coordinates

```bash
# Make GRCh37 gene annotation
magma --annotate \
  --snp-loc /Users/c24102394/MAGMA/ref/g1000_eur.bim \
  --gene-loc /Users/c24102394/MAGMA/ref/NCBI37.3.gene.loc \
  --out genes_b37
```


### 3. Re-calculate (per-trait and per-SNP pvals)

```bash
# Switch the files around for SCZ
python3 - <<'PY'
import gzip, math
from math import erf, sqrt

fin = gzip.open('/Users/c24102394/Desktop/PhD/AD_SZ_genes/LAVA_LDSC_nature/Data/AD.sumstats.gz','rt')
fout = open('AD_for_magma.txt','w')

header = next(fin)  # SNP  N  Z  A1  A2
fout.write('SNP P N\n')

for line in fin:
    snp, N, Z, A1, A2 = line.split()
    z = float(Z)
    Phi = 0.5*(1.0 + erf(abs(z)/sqrt(2.0)))
    p = 2.0*(1.0 - Phi)
    fout.write(f"{snp} {p:.6g} {N}\n")

fout.close()
PY
```


### 4. Run gene analysis for trait 1 (Alzheimer´s Disease)

```bash
magma \
  --bfile /Users/c24102394/MAGMA/ref/g1000_eur \
  --pval AD_for_magma.txt use=SNP,P ncol=N \
  --gene-annot genes_b37.genes.annot \
  --out AD_magma
```

### 5. Run gene analysis for trait 2 (Schizophrenia)

```bash
magma \
  --bfile /Users/c24102394/MAGMA/ref/g1000_eur \
  --pval SZ_for_magma.txt use=SNP,P ncol=N \
  --gene-annot genes_b37.genes.annot \
  --out SZ_magma
```

```bash
head AD_magma.genes.out
head SZ_magma.genes.out
```

```bash
# GENE       CHR      START       STOP  NSNPS  NPARAM      N        ZSTAT            P
# 148398       1     859993     879961      3       1  58546      0.54806      0.29183
# 26155        1     879583     894679      8       3  58546      0.48183      0.31496
# 84069        1     901872     910488      1       1  58546      0.44129       0.3295
# 84808        1     910579     917473      1       1  58546    -0.034348       0.5137
# 9636         1     948847     949920      1       1  58546      -1.8145       0.9652
# 375790       1     955503     991499      9       3  58546      0.44975      0.32644
# 54991        1    1017198    1051736     14       4  58546      0.22966      0.40918
# 254173       1    1109286    1133315     13       4  58546       1.8364     0.033153
# 8784         1    1138888    1142163      4       2  58546      0.81201      0.20839
```


### 6. Map gene number to gene NAME

```bash
python3 - <<'PY'
import pandas as pd
from pathlib import Path

GENES = Path("/Users/c24102394/MAGMA/ref") # NCBI37.3.gene.loc
OUT = Path("/Users/c24102394/Desktop/PhD/AD_SZ_genes/LAVA_LDSC_nature/outputs/lava/magma_runs") # AD_magma.genes.out / SZ_magma.genes.out

gene_loc_df = pd.read_csv(GENES / "NCBI37.3.gene.loc",
                          sep=r"\s+",
                          header=None)

ad_df = pd.read_csv(OUT / "AD_magma.genes.out",
                    sep=r"\s+",
                    header=0)

sz_df = pd.read_csv(OUT / "SZ_magma.genes.out",
                    sep=r"\s+",
                    header=0)

gene_loc_df.rename(columns={
    0: "GENE",
    1: "CHR",
    2: "START",
    3: "END",
    4: "STRAND",
    5: "GENE_ID"
}, inplace=True)

def add_gene_id(magma_df: pd.DataFrame, gene_loc_df: pd.DataFrame) -> pd.DataFrame:
    magma_df = magma_df.copy()
    gene_loc_df = gene_loc_df.copy()
    magma_df["GENE"] = magma_df["GENE"].astype(int)
    gene_loc_df["GENE"] = gene_loc_df["GENE"].astype(int)
    merged = magma_df.merge(
        gene_loc_df[["GENE", "GENE_ID"]],
        on="GENE",
        how="left"
    )
    return merged

if __name__ == "__main__":
    ad_df = add_gene_id(ad_df, gene_loc_df)
    sz_df = add_gene_id(sz_df, gene_loc_df)
    ad_df.to_csv(OUT / "AD_magma.genes.mapped.tsv", index=False, sep="\t")
    sz_df.to_csv(OUT / "SZ_magma.genes.mapped.tsv", index=False, sep="\t")
PY
```

```bash
head AD_magma.genes.mapped.tsv
head SZ_magma.genes.mapped.tsv
```

```bash
# GENE	CHR	START	STOP	NSNPS	NPARAM	N	ZSTAT	P	GENE_ID
# 148398	1	859993	879961	3	1	58546	0.54806	0.29183	SAMD11
# 26155	1	879583	894679	8	3	58546	0.48183	0.31496	NOC2L
# 84069	1	901872	910488	1	1	58546	0.44129	0.3295	PLEKHN1
# 84808	1	910579	917473	1	1	58546	-0.0343	0.5137	PERM1
# 9636	1	948847	949920	1	1	58546	-1.8145	0.9652	ISG15
# 375790	1	955503	991499	9	3	58546	0.44975	0.32644	AGRN
# 54991	1	1017198	1051736	14	4	58546	0.22966	0.40918	C1orf159
# 254173	1	1109286	1133315	13	4	58546	1.8364	0.03315	TTLL10
# 8784	1	1138888	1142163	4	2	58546	0.81201	0.20839	TNFRSF18
```


### 7. Reformat for genes of interest per phenotype

```bash
python3 - <<'PY'
import pandas as pd
from pathlib import Path

OUT = Path("/Users/c24102394/Desktop/PhD/AD_SZ_genes/LAVA_LDSC_nature/outputs/lava/magma_runs")

def get_exact_genes(lava_genes: pd.DataFrame, phenotype_genes: pd.DataFrame) -> pd.DataFrame:
    phenotype_genes = phenotype_genes.copy()
    lava_genes = lava_genes.copy()
    lava_genes.loc[:, lava_genes.columns.isin(["symbol", "CHR", "start", "stop"])]
    lava_genes.rename(columns={
        "symbol": "GENE_ID",
        "start": "START",
        "stop": "STOP"},
        inplace=True
    )
    merged = phenotype_genes.merge(
        lava_genes,
        on="GENE_ID",
        how="inner"
    )
    merged.drop(
        columns={
            "CHR_y",
            "START_y",
            "STOP_y"
        },
        inplace=True
    )
    return merged

def reformat_final_df(phenotype1: pd.DataFrame, phenotype2: pd.DataFrame, name1: str = "AD", name2: str = "SZ") -> pd.DataFrame:
    phenotype1 = phenotype1.copy()
    phenotype2 = phenotype2.copy()
    cols_to_rename1 = {col: f"{col}_{name1}" for col in phenotype1.columns if col != "GENE_ID"}
    cols_to_rename2 = {col: f"{col}_{name2}" for col in phenotype2.columns if col != "GENE_ID"}
    phenotype1.rename(columns=cols_to_rename1, inplace=True)
    phenotype2.rename(columns=cols_to_rename2, inplace=True)
    merged = phenotype1.merge(phenotype2, on="GENE_ID", how="inner")
    return merged

if __name__ == "__main__":
    genes_lava = pd.read_csv(
        "/Users/c24102394/ref/magma/MAPPINGS/lava_gene_summary_with_names.csv",
        sep=",",
        index_col=0
    )

    genes_magma_ad = pd.read_csv(
        "/Users/c24102394/Desktop/PhD/AD_SZ_genes/LAVA_LDSC_nature/outputs/lava/magma_runs/AD_magma.genes.mapped.tsv",
        sep="\t",
        index_col=0
    )

    genes_magma_sz = pd.read_csv(
        "/Users/c24102394/Desktop/PhD/AD_SZ_genes/LAVA_LDSC_nature/outputs/lava/magma_runs/SZ_magma.genes.mapped.tsv",
        sep="\t",
        index_col=0
    )

    ad = get_exact_genes(genes_lava, genes_magma_ad)
    sz = get_exact_genes(genes_lava, genes_magma_sz)
    ad.to_csv(OUT / "AD_final_mapping.tsv", sep="\t", index=False)
    sz.to_csv(OUT / "SZ_final_mapping.tsv", sep="\t", index=False)
    combined = reformat_final_df(ad, sz, name1="AD", name2="SZ")
    combined.to_csv(OUT / "AD_SZ_combined.tsv", sep="\t", index=False)
    print("Combined AD+SZ mapping saved to:", OUT / "AD_SZ_combined.tsv")
PY
```




