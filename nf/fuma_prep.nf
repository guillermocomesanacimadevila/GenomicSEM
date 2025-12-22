#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.out_base = "${workflow.launchDir}/outputs/fuma_prep"
params.log_base = "${workflow.launchDir}/logs/fuma_prep"
params.clump_base = "${workflow.launchDir}/outputs/clumping"
params.post_base = "${workflow.launchDir}/outputs/fuma_post"

process PREP_4_FUMA {
  // mkdir -p logs/fuma_prep
  // mkdir -p outputs/fuma_prep
  // mkdir -p outputs/fuma_prep/AD-SCZ
  // mkdir -p outputs/fuma_prep/AD-SCZ-LON
  // copy from clumping outputs each respective locus for each respetcie pheno
  // mkdir -p output/fuma_results
  // gzip the locus 0/1 file per trait
  
  tag "${set_id}"

  publishDir "${params.out_base}", mode: 'copy'

  input:
  tuple val(set_id), val(clump_dir), val(traits)

  output:
  path("${set_id}")

  script:
  def traits_str = traits.join(' ')
  """
  set -e
  mkdir -p ${set_id}/loci
  mkdir -p ${set_id}/meta
  echo -e "set_id\ttrait\tlocus_id\tfile_gz" > ${set_id}/meta/manifest.tsv

  for t in ${traits_str}; do
    mkdir -p ${set_id}/loci/\$t
    for f in ${clump_dir}/locus_*_\$t.tsv; do
      [ -f "\$f" ] || continue
      bn=\$(basename "\$f")
      locus_id=\$(echo "\$bn" | cut -d'_' -f2)
      out=${set_id}/loci/\$t/\${bn}.gz

      python3 - "\$f" "\$out" <<'PY'
import sys, gzip

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as f, gzip.open(outfile, 'wt') as out:
    header = f.readline().rstrip('\\n').split('\\t')
    out.write('\\t'.join(header) + '\\n')

    idx = {name: i for i, name in enumerate(header)}

    for line in f:
        row = line.rstrip('\\n').split('\\t')
        for col in ('N', 'CHR', 'POS'):
            if col in idx and row[idx[col]] != '':
                row[idx[col]] = str(int(float(row[idx[col]])))
        out.write('\\t'.join(row) + '\\n')
PY

      echo -e "${set_id}\t\$t\t\$locus_id\t\$out" >> ${set_id}/meta/manifest.tsv
    done
  done

  [ -f ${clump_dir}/lead_snps.tsv ] && cp ${clump_dir}/lead_snps.tsv ${set_id}/meta/
  [ -f ${clump_dir}/locus_coords.tsv ] && cp ${clump_dir}/locus_coords.tsv ${set_id}/meta/
  """
}

process FUMA_POST_DIRS {

  tag "${trait}"
  publishDir "${params.post_base}", mode: 'copy'

  input:
  tuple val(trait), val(clump_dir)

  output:
  path("${trait}")

  script:
  """
  set -e
  mkdir -p ${trait}
  for f in ${clump_dir}/locus_*_${trait}.tsv; do
    [ -f "\$f" ] || continue
    bn=\$(basename "\$f")
    locus_id=\$(echo "\$bn" | cut -d'_' -f2)
    mkdir -p ${trait}/locus_\${locus_id}
  done
  """
}

workflow {

  Channel.of(
    tuple("AD-SCZ", "${params.clump_base}/AD_SCZ", ["AD","SCZ"]),
    tuple("AD-SCZ-LON", "${params.clump_base}/AD-SCZ-LON", ["AD","SCZ","LON"])
  ) | PREP_4_FUMA | view

  Channel.of(
    tuple("AD", "${params.clump_base}/AD_SCZ"),
    tuple("SCZ", "${params.clump_base}/AD_SCZ"),
    tuple("AD", "${params.clump_base}/AD-SCZ-LON"),
    tuple("SCZ", "${params.clump_base}/AD-SCZ-LON"),
    tuple("LON", "${params.clump_base}/AD-SCZ-LON")
  ) | FUMA_POST_DIRS | view
}
