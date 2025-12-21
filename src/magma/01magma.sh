#!/usr/bin/env bash
set -euo pipefail

command -v magma >/dev/null 2>&1 || {
  echo "ERROR: MAGMA is not installed"
  exit 1
}
# check if magma insalled - else print("MAGMA != installed")

SNP_LOC=$1
GENE_LOC=$2
OUT_FILE=$3
OUT_DIR=$4

# magma --annotate \
  #  --snp-loc /Users/c24102394/MAGMA/ref/g1000_eur.bim \
  #  --gene-loc /Users/c24102394/MAGMA/ref/NCBI37.3.gene.loc \
  #  --out genes_b37

mkdir -p "$OUT"

magma --annotate \
  --snp-loc "$SNP_LOC" \
  --gene-loc "$GENE_LOC" \
  --out "$OUT_DIR/$OUT_FILE"

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



# Gene analysis for trait1:

# gene analysis for trait2

# map gene ID to gene name from GRCh37 on both

# 2 MAGMA mapped sumstats (1 per trait)

