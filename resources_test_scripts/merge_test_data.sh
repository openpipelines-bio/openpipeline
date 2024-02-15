
#!/bin/bash



# settings
ID=merge_test_data
OUT=resources_test/$ID
DIR="$OUT"

mkdir -p "$OUT"

target/docker/dataflow/split_modalities/split_modalities \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu \
  --output "$OUT"

