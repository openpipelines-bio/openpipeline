#!/bin/bash

set -eo pipefail

## VIASH START
par_input='https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5'
par_output='pbmc_1k_protein_v3_raw_feature_bc_matrix.h5'
par_verbose='false'
## VIASH END

extra_params=()

if [ "$par_verbose" != "true" ]; then
  extra_params+=("--quiet")
fi

wget "$par_input" -O "$par_output" "${extra_params[@]}"
