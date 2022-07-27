
#!/bin/bash

# settings
ID=merge
OUT=resources_test/$ID
DIR="$OUT"

mkdir -p "$OUT"

target/docker/split/split_modalities/split_modalities \
--input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu \
--output "$OUT"

