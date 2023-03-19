#!/bin/bash
set -eo pipefail

# add additional params
extra_params=( )

if [ ! -z "$par_sm" ]; then
  extra_params+=( "--sm $par_sm" )
fi

if [ ! -z "$par_smList" ]; then 
  extra_params+=( "--sm-list $par_smList" )
fi

if [ $par_skipUmi == "true" ]; then 
  extra_params+=( "--skip-umi" )
fi

if [ ! -z "$par_groupList" ]; then 
  extra_params+=( "--group-list $par_groupList" )
fi

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi


popscle dsc-pileup --sam $par_sam --tag-group $par_tagGroup --tag-UMI $par_tagUMI \
       --exclude-flag $par_excludeFlag --sam-verbose $par_samVerbose \
       --vcf $par_vcf --vcf-verbose $par_vcfVerbose --out ${par_output}/${par_out} \
       --cap-BQ $par_capBQ --min-BQ $par_minBQ --min-MQ $par_minMQ --min-TD $par_minTD \
       --excl-flag $par_exclFlag --min-total $par_minTotal --min-uniq $par_minUniq \
       --min-snp $par_minSnp ${extra_params[@]}
