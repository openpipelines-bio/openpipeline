#!/bin/bash

set -eo pipefail

# Unset flags if they equal 'false'
[[ "$par_skip_remap" == "false" ]] && unset par_skip_remap
[[ "$par_ignore" == "false" ]] && unset par_ignore
[[ "$par_souporcell3" == "false" ]] && unset par_souporcell3
[[ "$par_no_umi" == "false" ]] && unset par_no_umi

# souporcell only creates the logs directory when the output directory
# does not exist yet, so make sure it is present.
mkdir -p "$par_output/logs"

known_genotypes_sample_names=()
if [ -n "$par_known_genotypes_sample_names" ]; then
  IFS=';' read -r -a known_genotypes_sample_names <<< "$par_known_genotypes_sample_names"
fi

# Note: souporcell parses boolean options with python's bool(), so any
# non-empty value (including 'False') enables them. Only pass them when set.
souporcell_pipeline.py \
  --bam "$par_bam" \
  --fasta "$par_fasta" \
  --barcodes "$par_barcodes" \
  --clusters "$par_clusters" \
  --ploidy "$par_ploidy" \
  --min_alt "$par_min_alt" \
  --min_ref "$par_min_ref" \
  --max_loci "$par_max_loci" \
  --umi_tag "$par_umi_tag" \
  --out_dir "$par_output" \
  --threads "${meta_cpus:-1}" \
  ${par_clustering_method:+--clustering_method "$par_clustering_method"} \
  ${par_souporcell3:+--souporcell3 True} \
  ${par_restarts:+--restarts "$par_restarts"} \
  ${par_common_variants:+--common_variants "$par_common_variants"} \
  ${par_known_genotypes:+--known_genotypes "$par_known_genotypes"} \
  ${par_known_genotypes_sample_names:+--known_genotypes_sample_names "${known_genotypes_sample_names[@]}"} \
  ${par_max_base_mem:+--max_base_mem "$par_max_base_mem"} \
  ${par_no_umi:+--no_umi True} \
  ${par_skip_remap:+--skip_remap True} \
  ${par_ignore:+--ignore True}

cut -d$'\t' -f 1-3 "$par_output/clusters.tsv" |
sed 's/\t/,/g' |
awk 'BEGIN{FS=OFS=","} {$2=($2=="singlet")?$3:$2; NF=NF-1; print}' |
sed '1s/barcode,status/cell,donor_id/' > "$par_output/cell_annotation.csv"
