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

