#!/bin/bash

# Creates additional h5mu files for testing the components based on pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu
# You need to run src/download/sync_test_resources first to download the inital files.

# Get the root of the directory and enter it
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# Path variables
IN=resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu
OUT_FOLDER=resources_test/pbmc_1k_protein_v3/
OUT_NEIGHBORS="${OUT_FOLDER}pbmc_1k_protein_v3_filtered_feature_bc_matrix.neighbors.h5mu"
OUT_NEIGHBORS_LEIDEN="${OUT_FOLDER}pbmc_1k_protein_v3_filtered_feature_bc_matrix.neighbors.leiden.h5mu"

echo "Performing neighbors on ${IN}"

bin/viash run src/neighbors/find_neighbors/config.vsh.yaml -- \
    --input "${IN}" \
    --output "${OUT_NEIGHBORS}"

echo "Performing leiden on ${OUT_NEIGHBORS}"

bin/viash run src/cluster/leiden/config.vsh.yaml -- \
    --input "${OUT_NEIGHBORS}" \
    --output "${OUT_NEIGHBORS_LEIDEN}"
