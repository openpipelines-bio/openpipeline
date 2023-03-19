#!/bin/bash

set -eo pipefail

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

popscle freemuxlet \
  --plp $par_plp \
  --nsample $par_nsample \
  --verbose $par_verbose \
  --doublet-prior $par_doubletPrior \
  --geno-error $par_genoError \
  --bf-thres $par_bfThres \
  --frac-init-clust $par_fracInitClust \
  --iter-init $par_iterInit \
  --seed $par_seed \
  --cap-BQ $par_capBQ \
  --min-BQ $par_minBQ \
  --min-total $par_minTotal \
  --min-umi $par_minUmi \
  --min-snp $par_minSnp \
  --out "$par_output/${par_out}" \
  ${par_initCluster:+--init-cluster $par_initCluster} \
  ${par_auxFiles:+--aux-files} \
  ${par_keepInitMissing:+--keep-init-missing} \
  ${par_randomizeSingletScore:+--randomize-singlet-score} \
  ${par_groupList:+--group-list $par_groupList}

Rscript summary.R --freemuxlet_out "${par_output}/${par_out}.clust1.samples.gz"
