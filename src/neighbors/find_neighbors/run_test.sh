#!/usr/bin/env bash

set -ex

./find_neighbors --input pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad --output output.h5ad
[[ ! -f output.h5ad ]] && echo "Output file could not be found!" && exit 1

echo ">>> Test finished successfully"
