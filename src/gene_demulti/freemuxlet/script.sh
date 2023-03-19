#!/bin/bash

set -eo pipefail

# add additional params
extra_params=( )

if [ ! -z "$par_initCluster" ]; then 
  extra_params+=( "--init-cluster $par_initCluster" )
fi

if [ $par_auxFiles == true ]; then
  extra_params+=( "--aux-files" )
fi

if [ $par_keepInitMissing == true ]; then
  extra_params+=( "--keep-init-missing" )
fi

if [ $par_randomizeSingletScore == true ]; then
  extra_params+=( "--randomize-singlet-score" )
fi

if [ ! -z "$par_groupList" ]; then
  extra_params+=( "--group-list $par_groupList" )
fi

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

popscle freemuxlet --plp $par_plp --nsample $par_nsample --verbose $par_verbose \
       --doublet-prior $par_doubletPrior --geno-error $par_genoError --bf-thres $par_bfThres \
       --frac-init-clust $par_fracInitClust --iter-init $par_iterInit  \
       --seed $par_seed --cap-BQ $par_capBQ --min-BQ $par_minBQ \
       --min-total $par_minTotal --min-umi $par_minUmi --min-snp $par_minSnp \
       --out $par_output${par_out} ${extra_params[@]}

Rscript summary.R --freemuxlet_out ${par_output}/${par_out}.clust1.samples.gz
