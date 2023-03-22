#!/bin/bash
set -eo pipefail

# Unset flags if they equal 'false'
[[ "$par_skipUmi" == "false" ]] && unset par_skipUmi

# Create output directory if it doesn't exist
if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

popscle dsc-pileup \
  --sam $par_sam \
  --tag-group $par_tagGroup \
  --tag-UMI $par_tagUMI \
  --exclude-flag $par_excludeFlag \
  --sam-verbose $par_samVerbose \
  --vcf $par_vcf \
  --vcf-verbose $par_vcfVerbose \
  --out "$par_output/$par_out" \
  --cap-BQ $par_capBQ \
  --min-BQ $par_minBQ \
  --min-MQ $par_minMQ \
  --min-TD $par_minTD \
  --excl-flag $par_exclFlag \
  --min-total $par_minTotal \
  --min-uniq $par_minUniq \
  --min-snp $par_minSnp \
  ${par_sm:+--sm $par_sm} \
  ${par_smList:+--sm-list $par_smList} \
  ${par_skipUmi:+--skip-umi} \
  ${par_groupList:+--group-list $par_groupList}
