#!/bin/bash
set -eo pipefail

# Unset flags if they equal 'false'
[[ "$par_skip_umi" == "false" ]] && unset par_skip_umi

# Create output directory if it doesn't exist
if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

popscle dsc-pileup \
  --sam $par_sam \
  --tag-group $par_tag_group \
  --tag-UMI $par_tag_umi \
  --exclude-flag $par_exclude_flag \
  --sam-verbose $par_sam_verbose \
  --vcf $par_vcf \
  --vcf-verbose $par_vcf_verbose \
  --out "$par_output/$par_out" \
  --cap-BQ $par_cap_bq \
  --min-BQ $par_min_bq \
  --min-MQ $par_min_mq \
  --min-TD $par_min_td \
  --excl-flag $par_excl_flag \
  --min-total $par_min_total \
  --min-uniq $par_min_uniq \
  --min-snp $par_min_snp \
  ${par_sm:+--sm $par_sm} \
  ${par_sm_list:+--sm-list $par_sm_list} \
  ${par_skip_umi:+--skip-umi} \
  ${par_group_list:+--group-list $par_group_list}
