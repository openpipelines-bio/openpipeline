#!/usr/bin/env bash

set -ex

./lognorm --input pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu --output output.h5mu
[[ ! -f output.h5mu ]] && echo "Output file could not be found!" && exit 1

echo ">>> Test finished successfully"
