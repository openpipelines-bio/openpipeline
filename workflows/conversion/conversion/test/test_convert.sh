#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

# No input_type, should fail
nextflow \
  run . \
  -main-script workflows/1_ingestion/conversion/main.nf \
  --input ./resources_test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  --output output/ \
  -resume \
  -c ./workflows/1_ingestion/conversion/nextflow.config

# Standard 10x h5 conversion
nextflow \
  run . \
  -main-script workflows/1_ingestion/conversion/main.nf \
  --input_type 10xh5 \
  --input ./resources_test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  --output output/ \
  -resume \
  -c ./workflows/1_ingestion/conversion/nextflow.config
