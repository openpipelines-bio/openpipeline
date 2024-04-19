#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --vcf "$meta_resources_dir/demuxafy_test_data/test_dataset_chr1_2.vcf" \
    --vcf "$meta_resources_dir/demuxafy_test_data/test_dataset_chr3_4.vcf" \
    --concat --filter \
    --output bcftools_result/

[[ ! -f bcftools_result/filtered_sorted_concated_chroms.vcf ]] && echo "Output processed VCF file could not be found!" && exit 1

echo ">>> Test finished successfully"
