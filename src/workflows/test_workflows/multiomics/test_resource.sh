#!/bin/bash

set -e

viash run src/dataflow/split_modalities/config.vsh.yaml -- \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu \
  --output output_test/split_modalities/h5mu/ \
  --output_types output_test/split_modalities/foo_types.csv \
  --compression gzip