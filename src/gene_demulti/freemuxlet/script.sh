#!/bin/bash

set -eo pipefail

# add additional params
extra_paramsDsc=( )
extra_paramsFreemuxlet=( )
extra_params=( )

if [ ! -z "$par_plp" ] && [ $par_plp != "true" ]; then 
  extra_paramsFreemuxlet+=( "--plp $par_plp" )
fi
  
if [ ! -z "$par_sm" ]; then 
  extra_paramsDsc+=( "--sm $par_sm" )
fi

if [ ! -z "$par_smList" ]; then 
  extra_paramsDsc+=( "--sm-list $par_smList" )
fi

if [ $par_skipUmi == true ]; then
  extra_paramsDsc+=( "--skip-umi" )
fi

if [ ! -z "$par_groupList" ]; then 
  extra_params+=( "--group-list $par_groupList" )
fi

if [ ! -z "$par_initCluster" ]; then 
  extra_paramsFreemuxlet+=( "--init-cluster $par_initCluster" )
fi

if [ $par_auxFiles == true ]; then
  extra_paramsFreemuxlet+=( "--aux-files" )
fi

if [ $par_keepInitMissing == true ]; then
  extra_paramsFreemuxlet+=( "keep-init-missing" )
fi

if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi

if [! -z "$par_plp" ] && [ "$par_plp" != "true" ]; then
  popscle freemuxlet --verbose $par_verbose --doublet-prior $par_doubletPrior \
       --bf-thres $par_bfThres --cap-BQ $par_capBQ --min-BQ $par_minBQ  \
       --frac-init-clust $par_fracInitClust --iter-init $par_iterInit \
       --min-total $par_minTotal --min-umi $par_minUmi --min-snp $par_minSnp \
       --out ${par_output}${par_out} ${extra_params[@]} ${extra_paramsFreemuxlet[@]}
else
  popscle dsc-pileup --sam $par_sam --tag-group $par_tagGroup --tag-UMI $par_tagUMI \
       --exclude-flag $par_excludeFlag --sam-verbose $par_samVerbose \
       --vcf $par_vcf --vcf-verbose $par_vcfVerbose --out ${par_output}${par_outDsc} \
       --cap-BQ $par_capBQDsc --min-BQ $par_minBQ --min-MQ $par_minMQ --min-TD $par_minTD \
       --excl-flag $par_exclFlag --min-total $par_minTotal --min-uniq $par_minUniq \
       --min-snp $par_minSnp ${extra_params[@]} ${extra_paramsDsc[@]} 

  popscle freemuxlet --frac-init-clust $par_fracInitClust --iter-init $par_iterInit \
       --verbose $par_verbose --doublet-prior $par_doubletPrior --bf-thres $par_bfThres \
       --cap-BQ $par_capBQ --min-BQ $par_minBQ  --plp ${par_output}${par_outDsc} \
       --min-total $par_minTotal --min-umi $par_minUmi --min-snp $par_minSnp \
       --out $par_output${par_out} ${extra_params[@]} ${extra_paramsFreemuxlet[@]}
fi
