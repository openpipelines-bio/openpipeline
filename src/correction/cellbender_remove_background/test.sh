
#!/bin/bash

echo "> Running CellBender"

$meta_executable \
  --input "$meta_resources_dir/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu" \
  --output "foo.h5mu" \
  --epochs 100

[ -f foo.h5mu ] || { echo "No foo.h5mu output file found!"; exit 1; }

# todo: upgrade test script to python, read in h5mu and test contents


echo "> Test succeeded!"