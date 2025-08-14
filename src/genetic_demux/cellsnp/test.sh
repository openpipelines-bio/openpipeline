#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --sam_file "$meta_resources_dir/demuxafy_test_data/chr_1_pooled.sorted.bam" \
    --barcode_file "$meta_resources_dir/demuxafy_test_data/barcodes.tsv" \
    --output cellsnp_result/ \
    --regions_vcf "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf"

[[ ! -f cellsnp_result/cellSNP.base.vcf ]] && echo "Output VCF file could not be found!" && exit 1

echo ">>> Test finished successfully"
