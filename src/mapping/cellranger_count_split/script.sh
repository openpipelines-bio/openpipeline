#!/bin/bash

set -eo pipefail

## VIASH START
par_input='input_dir'
par_filtered_h5='filtered_feature_bc_matrix.h5'
par_metrics_summary='metrics_summary.csv'
par_molecule_info='molecule_info.h5'
par_bam='possorted_genome_bam.bam'
par_bai='possorted_genome_bam.bam.bai'
par_raw_h5='raw_feature_bc_matrix.h5'
## VIASH END

filtered_h5="$par_input/filtered_feature_bc_matrix.h5"
if [ -f "$filtered_h5" ] && [ ! -z "$par_filtered_h5" ]; then
  echo "+ cp $filtered_h5 $par_filtered_h5"
  cp "$filtered_h5" "$par_filtered_h5"
fi

metrics_summary="$par_input/metrics_summary.csv"
if [ -f "$metrics_summary" ] && [ ! -z "$par_metrics_summary" ]; then
  echo "+ cp $metrics_summary $par_metrics_summary"
  cp "$metrics_summary" "$par_metrics_summary"
fi

molecule_info="$par_input/molecule_info.h5"
if [ -f "$molecule_info" ] && [ ! -z "$par_molecule_info" ]; then
  echo "+ cp $molecule_info $par_molecule_info"
  cp "$molecule_info" "$par_molecule_info"
fi

bam="$par_input/possorted_genome_bam.bam"
if [ -f "$bam" ] && [ ! -z "$par_bam" ]; then
  echo "cp $bam $par_bam"
  cp "$bam" "$par_bam"
fi

raw_h5="$par_input/raw_feature_bc_matrix.h5"
if [ -f "$raw_h5" ] && [ ! -z "$par_raw_h5" ]; then
  echo "+ cp $raw_h5 $par_raw_h5"
  cp "$raw_h5" "$par_raw_h5"
fi

bai="$par_input/possorted_genome_bam.bam.bai"
if [ -f "$bai" ] && [ ! -z "$par_bai" ]; then
  echo "+ cp $bai $par_bai"
  cp "$bai" "$par_bai"
fi