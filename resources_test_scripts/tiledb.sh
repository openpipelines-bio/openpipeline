#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=tiledb
OUT="resources_test/$ID"

# create raw directory
raw_dir="$OUT/"
mkdir -p "$raw_dir"


viash run src/convert/from_h5mu_or_h5ad_to_tiledb/config.vsh.yaml -- \
    --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu \
    --tiledb_dir "$OUT"/pbmc_1k_protein_v3_mms \
    --rna_modality "rna" \
    --rna_raw_layer_input "X" \
    --rna_normalized_layer_input "log_normalized" \
    --rna_var_gene_names_input "gene_symbol"