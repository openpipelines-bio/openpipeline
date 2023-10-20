#!/bin/bash

set -eo pipefail

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

popscle freemuxlet \
  --plp $par_plp \
  --nsample $par_nsample \
  --verbose $par_verbose \
  --doublet-prior $par_doublet_prior \
  --geno-error $par_geno_error \
  --bf-thres $par_bf_thres \
  --frac-init-clust $par_frac_init_clust \
  --iter-init $par_iter_init \
  --seed $par_seed \
  --cap-BQ $par_cap_bq \
  --min-BQ $par_min_bq \
  --min-total $par_min_total \
  --min-umi $par_min_umi \
  --min-snp $par_min_snp \
  --out "$par_output/$par_out" \
  ${par_init_cluster:+--init-cluster $par_init_cluster} \
  ${par_aux_files:+--aux-files} \
  ${par_keep_init_missing:+--keep-init-missing} \
  ${par_randomize_singlet_score:+--randomize-singlet-score} \
  ${par_group_list:+--group-list $par_group_list}

Rscript summary.R --freemuxlet_out "$par_output/$par_out.clust1.samples.gz"
