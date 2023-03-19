#!/bin/bash
set -eo pipefail

# Create output directory if it doesn't exist
if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

popscle demuxlet \
  --tag-group $par_tagGroup \
  --tag-UMI $par_tagUMI \
  --field $par_field \
  --geno-error-offset $par_genoErrorOffset \
  --geno-error-coeff $par_genoErrorCoeff \
  --r2-info $par_r2Info \
  --min-mac $par_minMac \
  --min-callrate $par_minCallrate \
  --alpha $par_alpha \
  --doublet-prior $par_doubletPrior \
  --sam-verbose $par_samVerbose \
  --vcf $par_vcf \
  --vcf-verbose $par_vcfVerbose \
  --cap-BQ $par_capBQ \
  --min-BQ $par_minBQ \
  --min-MQ $par_minMQ \
  --min-TD $par_minTD \
  --excl-flag $par_exclFlag \
  --min-total $par_minTotal \
  --min-umi $par_minUmi \
  --min-snp $par_minSnp \
  --out ${par_output}/${par_out} \
  ${par_plp:+--plp $par_plp} \
  ${par_sam:+--sam $par_sam} \
  ${par_sm:+--sm $par_sm} \
  ${par_smList:+--sm-list $par_smList} \
  ${par_groupList:+--group-list $par_groupList}

Rscript summary.R --demuxlet_out ${par_output}/${par_out}.best