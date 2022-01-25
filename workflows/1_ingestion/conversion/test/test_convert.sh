#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

# No input_type, should fail
nextflow \
  run . \
  -main-script workflows/1_ingestion/conversion/main.nf \
  --input ./resources/test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  --output output/ \
  -resume \
  -c ./workflows/1_ingestion/conversion/nextflow.config

# Standard 10x h5 conversion
nextflow \
  run . \
  -main-script workflows/1_ingestion/conversion/main.nf \
  --input_type 10x_h5 \
  --input ./resources/test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  --output output/ \
  -resume \
  -c ./workflows/1_ingestion/conversion/nextflow.config

# 10x mtx conversion
nextflow \
  run . \
  -main-script workflows/1_ingestion/conversion/main.nf \
  --input_type 10x_mtx \
  --input ./resources/test/convert/convert_10x_mtx_to_h5ad/ \
  --output output/ \
  -resume \
  -c ./workflows/1_ingestion/conversion/nextflow.config

# 10x csv conversion
nextflow \
  run . \
  -main-script workflows/1_ingestion/conversion/main.nf \
  --input_type csv \
  --input ./resources/test/convert/convert_csv_to_h5ad/CS0000007_subsample_LI00080.csv.gz \
  --output output/ \
  -resume \
  -c ./workflows/1_ingestion/conversion/nextflow.config

