#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --cell_data "$meta_resources_dir/vireo_test_data/cells.cellSNP.vcf.gz" \
    --n_donor 4 \
    --output vireo_result/

[[ ! -f vireo_result/GT_donors.vireo.vcf.gz ]] && echo "Output donor genotype file could not be found!" && exit 1
[[ ! -f vireo_result/cell_annotation.csv ]] && echo "Output cell annotation csv could not be found!" && exit 1
echo ">>> Test finished successfully"
