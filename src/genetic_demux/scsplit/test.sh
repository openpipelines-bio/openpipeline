#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --bam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --vcf "$meta_resources_dir/demuxafy_test_data/mixed_variant.vcf" \
    --num 12 \
    --ems 1 \
    --com "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf" \
    --bar "$meta_resources_dir/demuxafy_test_data/barcodes.tsv" \
    --ref ref_filtered.csv \
    --alt alt_filtered.csv \
    --output scSplit_result/

[[ ! -f scSplit_result/scSplit_result.csv ]] && echo "Output donor assignment file could not be found!" && exit 1
[[ ! -f scSplit_result/cell_annotation.csv ]] && echo "Output cell annotation file as tsv format could not be found!" && exit 1
echo ">>> Test finished successfully"
