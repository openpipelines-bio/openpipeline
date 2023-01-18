#!/bin/bash
set -eo pipefail

# add additional params
extra_params=( )
extra_paramsDemuxlet=( )
extra_paramsDsc=( )

if [ ! -z "$par_plp" ] && [ $par_plp != "true" ]; then 
  extra_paramsDemuxlet+=( "--plp $par_plp" )
fi

if [ -z "$par_plp" ]; then 
  extra_paramsDemuxlet+=( "--sam $par_sam" )
fi
  
if [ ! -z "$par_sm" ]; then 
  extra_params+=( "--sm $par_sm" )
fi

if [ ! -z "$par_smList" ]; then 
  extra_params+=( "--sm-list $par_smList" )
fi

if [ $par_skipUmi == "true" ]; then 
  extra_paramsDsc+=( "--skip-umi" )
fi

if [ ! -z "$par_groupList" ]; then 
  extra_params+=( "--group-list $par_groupList" )
fi

if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi


if [ -z "$par_plp" ] || [ "$par_plp" != "true" ]; then
  popscle demuxlet --tag-group $par_tagGroup --tag-UMI $par_tagUMI \
       --field $par_field --geno-error-offset $par_genoErrorOffset \
       --geno-error-coeff $par_genoErrorCoeff --r2-info $par_r2Info \
       --min-mac $par_minMac --min-callrate $par_minCallrate --alpha $par_alpha \
       --doublet-prior $par_doubletPrior --sam-verbose $par_samVerbose \
       --vcf $par_vcf --vcf-verbose $par_vcfVerbose --cap-BQ $par_capBQ --min-BQ $par_minBQ \
       --min-MQ $par_minMQ --min-TD $par_minTD --excl-flag $par_exclFlag \
       --min-total $par_minTotal --min-umi $par_minUmi --min-snp $par_minSnp \
       --out ${par_output}${par_out} ${extra_params[@]} ${extra_paramsDemuxlet[@]}

else
  popscle dsc-pileup --sam $par_sam --tag-group $par_tagGroup --tag-UMI $par_tagUMI \
       --exclude-flag $par_excludeFlag --sam-verbose $par_samVerbose \
       --vcf $par_vcfDsc --vcf-verbose $par_vcfVerbose --out ${par_output}${par_outDsc} \
       --cap-BQ $par_capBQDsc --min-BQ $par_minBQ --min-MQ $par_minMQ --min-TD $par_minTD \
       --excl-flag $par_exclFlag --min-total $par_minTotal --min-uniq $par_minUniq \
       --min-snp $par_minSnp ${extra_params[@]} ${extra_paramsDsc[@]} 

  popscle demuxlet --tag-group $par_tagGroup --tag-UMI $par_tagUMI --plp ${par_output}${par_outDsc} \
       --field $par_field --geno-error-offset $par_genoErrorOffset \
       --geno-error-coeff $par_genoErrorCoeff --r2-info $par_r2Info \
       --min-mac $par_minMac --min-callrate $par_minCallrate --alpha $par_alpha \
       --doublet-prior $par_doubletPrior --sam-verbose $par_samVerbose \
       --vcf $par_vcf --vcf-verbose $par_vcfVerbose --cap-BQ $par_capBQ --min-BQ $par_minBQ \
       --min-MQ $par_minMQ --min-TD $par_minTD --excl-flag $par_exclFlag \
       --min-total $par_minTotal --min-umi $par_minUmi --min-snp $par_minSnp \
       --out ${par_output}${par_out} ${extra_params[@]} ${extra_paramsDemuxlet[@]}
fi

Rscript summary.R --demuxlet_out ${par_output}${par_out}.best
