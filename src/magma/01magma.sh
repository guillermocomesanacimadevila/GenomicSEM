#!/usr/bin/env bash
set -euo pipefail

# chmod +x magma.sh && ./src/magma/magma.sh \
  #  /outputs/magma/AD/MAGMA \
  #  /Users/c24102394/MAGMA/ref/g1000_eur \
  #  /Users/c24102394/MAGMA/ref/g1000_eur.bim \
  #  /Users/c24102394/MAGMA/ref/NCBI37.3.gene.loc \
  #  AD \
  #  /data/Main/AD/post-ldsc/AD.sumstats.gz

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 OUT_DIR BFILE_PREFIX SNP_LOC GENE_LOC TRAIT_NAME SUMSTATS_GZ" >&2
  exit 1
fi

OUT_DIR="$1"
BFILE_PREFIX="$2"
SNP_LOC="$3"
GENE_LOC="$4"
TRAIT="$5"
SUMSTATS="$6"

MAGMA_BIN="${MAGMA_BIN:-/Users/c24102394/MAGMA/magma}"

if [[ ! -x "$MAGMA_BIN" ]]; then
  echo "ERROR: MAGMA binary not executable at: $MAGMA_BIN" >&2
  exit 1
fi

if [[ "$(uname -s)" == "Darwin" ]]; then
  MAGMA_CMD=(arch -x86_64 "$MAGMA_BIN")
else
  MAGMA_CMD=("$MAGMA_BIN")
fi

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

if [[ ! -f genes_b37.genes.annot ]]; then
  "${MAGMA_CMD[@]}" --annotate --snp-loc "$SNP_LOC" --gene-loc "$GENE_LOC" --out genes_b37
fi

PV_FILE="${TRAIT}_for_magma.txt"

python3 - --input "$SUMSTATS" --output "$PV_FILE" <<'PY'
import argparse
import gzip
from math import erf, sqrt

def z_to_p(z: float) -> float:
    Phi = 0.5 * (1.0 + erf(abs(z) / sqrt(2.0)))
    return 2.0 * (1.0 - Phi)

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()

def run(a: argparse.Namespace) -> None:
    opener = gzip.open if a.input.endswith(".gz") else open
    with opener(a.input, "rt") as fin, open(a.output, "w") as fout:
        next(fin, None)
        fout.write("SNP P N\n")
        for line in fin:
            if not line.strip():
                continue
            snp, n, z, *_ = line.split()
            fout.write(f"{snp} {z_to_p(float(z)):.6g} {n}\n")

def main() -> None:
    run(parse_args())

if __name__ == "__main__":
    main()
PY

"${MAGMA_CMD[@]}" --bfile "$BFILE_PREFIX" --pval "$PV_FILE" use=SNP,P ncol=N --gene-annot genes_b37.genes.annot --out "${TRAIT}_magma"

MAPPED_TSV="${TRAIT}_magma.genes.mapped.tsv"

python3 - --gene-loc "$GENE_LOC" --magma-out "${TRAIT}_magma.genes.out" --out-tsv "$MAPPED_TSV" <<'PY'
import argparse
import pandas as pd

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--gene-loc", required=True)
    p.add_argument("--magma-out", required=True)
    p.add_argument("--out-tsv", required=True)
    return p.parse_args()

def run(a: argparse.Namespace) -> None:
    gene_loc = pd.read_csv(a.gene_loc, sep=r"\s+", header=None).rename(columns={0: "GENE", 5: "GENE_ID"})
    gene_loc["GENE"] = gene_loc["GENE"].astype(int)
    magma = pd.read_csv(a.magma_out, sep=r"\s+", header=0)
    magma["GENE"] = magma["GENE"].astype(int)
    out = magma.merge(gene_loc[["GENE", "GENE_ID"]], on="GENE", how="left")
    out.to_csv(a.out_tsv, sep="\t", index=False)

def main() -> None:
    run(parse_args())

if __name__ == "__main__":
    main()
PY

echo "${TRAIT}_magma.genes.out"
echo "${TRAIT}_magma.genes.mapped.tsv"
