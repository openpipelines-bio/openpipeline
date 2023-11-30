#!/bin/bash

set -eo pipefail

# Unset flags if they equal 'false'
[[ "$par_no_doublet" == "false" ]] && unset par_no_doublet
[[ "$par_force_learn_gt" == "false" ]] && unset par_force_learn_gt
[[ "$par_ase_mode" == "false" ]] && unset par_ase_mode
[[ "$par_no_plot" == "false" ]] && unset par_no_plot
[[ "$par_call_ambient_rnas" == "false" ]] && unset par_call_ambient_rnas

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

vireo \
  --cellData $par_cell_data \
  --nDonor $par_n_donor \
  --genoTag $par_geno_tag \
  --nInit $par_n_init \
  --extraDonor $par_extra_donor \
  --out "${par_output}" \
  ${par_vartrix_data:+--vatrixData $par_vartrix_data} \
  ${par_donor_file:+--donorFile $par_donor_file} \
  ${par_no_doublet:+--noDoublet} \
  ${par_extra_donorMode:+--extraDonorMode $par_extra_donorMode} \
  ${par_force_learn_gt:+--forceLearnGT} \
  ${par_ase_mode:+--ASEmode} \
  ${par_no_plot:+--noPlot} \
  ${par_rand_seed:+--randSeed $par_rand_seed} \
  ${par_cell_range:+--cellRange $par_cell_range} \
  ${par_call_ambient_rnas:+--callAmbientRNAs} \
  ${meta_cpus:+--nproc $meta_cpus}

cut -d$'\t' -f 1-2 "$par_output/donor_ids.tsv" | tr '\t' ',' > "$par_output/cell_annotation.csv"
awk 'BEGIN{FS=OFS=","} NR>1{ gsub("donor", "", $2) } 1' "$par_output/cell_annotation.csv" > "$par_output/cell_annotation_temp.csv" && mv "$par_output/cell_annotation_temp.csv" "$par_output/cell_annotation.csv"
