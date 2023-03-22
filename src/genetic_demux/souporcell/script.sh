#!/bin/bash

set -eo pipefail

# Unset flags if they equal 'false'
[[ "$par_skip_remap" == "false" ]] && unset par_skip_remap
[[ "$par_ignore" == "false" ]] && unset par_ignore

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

/opt/souporcell/souporcell_pipeline.py \
  --bam $par_bam \
  --fasta $par_fasta \
  --barcodes $par_barcodes \
  --clusters $par_clusters \
  --ploidy $par_ploidy \
  --min_alt $par_min_alt \
  --min_ref $par_min_ref \
  --max_loci $par_max_loci \
  --out_dir $par_output \
  --threads ${par_threads:=1} \
  ${par_restarts:+--restarts $par_restarts} \
  ${par_common_variants:+--common_variants $par_common_variants} \
  ${par_known_genotypes:+--known_genotypes $par_known_genotypes} \
  ${par_known_genotypes_sample_names:+--known_genotypes_sample_names $par_known_genotypes_sample_names} \
  ${par_skip_remap:+--skip_remap True} \
  ${par_ignore:+--ignore True}

cut -d$'\t' -f 1-2 "$par_output/clusters.tsv" > "$par_output/assignment.tsv"
sed 1d "$par_output/assignment.tsv" > $par_output/res.tsv
{ echo cell$'\t'donor_id; cat "$par_output/res.tsv"; } > "$par_output/assignment.tsv"
rm ${par_output}res.tsv
