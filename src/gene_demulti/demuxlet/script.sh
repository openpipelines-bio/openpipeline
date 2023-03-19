#!/bin/bash
set -eo pipefail

# add additional params
extra_params=( )

if [ ! -z "$par_plp" ]; then
  extra_params+=( "--plp $par_plp" )
fi

if [ ! -z "$par_sam" ]; then
  extra_params+=( "--sam $par_sam" )
fi
  
if [ ! -z "$par_sm" ]; then 
  extra_params+=( "--sm $par_sm" )
fi

if [ ! -z "$par_smList" ]; then 
  extra_params+=( "--sm-list $par_smList" )
fi

if [ ! -z "$par_groupList" ]; then 
  extra_params+=( "--group-list $par_groupList" )
fi

if [ ! -d "$par_output" ]; then
  mkdir -p $par_output
fi

popscle demuxlet --tag-group $par_tagGroup --tag-UMI $par_tagUMI \
       --field $par_field --geno-error-offset $par_genoErrorOffset \
       --geno-error-coeff $par_genoErrorCoeff --r2-info $par_r2Info \
       --min-mac $par_minMac --min-callrate $par_minCallrate --alpha $par_alpha \
       --doublet-prior $par_doubletPrior --sam-verbose $par_samVerbose \
       --vcf $par_vcf --vcf-verbose $par_vcfVerbose --cap-BQ $par_capBQ --min-BQ $par_minBQ \
       --min-MQ $par_minMQ --min-TD $par_minTD --excl-flag $par_exclFlag \
       --min-total $par_minTotal --min-umi $par_minUmi --min-snp $par_minSnp \
       --out ${par_output}/${par_out} ${extra_params[@]}

Rscript summary.R --demuxlet_out ${par_output}/${par_out}.best
