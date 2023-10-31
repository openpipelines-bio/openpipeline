#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --bam "$meta_resources_dir/demuxafy_test_data/chr_1_pooled.sorted.bam" \
    --bam_index "$meta_resources_dir/demuxafy_test_data/chr_1_pooled.sorted.bam.bai" \
    --fasta "$meta_resources_dir/demuxafy_test_data/genome_chr1.fa" \
    --barcodes "$meta_resources_dir/demuxafy_test_data/barcodes.tsv" \
    --common_variants "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf" \
    --clusters 14 \
    --skip_remap \
    --output soup_result/

[[ ! -f soup_result/clusters.tsv ]] && echo "Output donor assignment file could not be found!" && exit 1
[[ ! -f soup_result/cell_annotation.csv ]] && echo "Output cell annotation file could not be found!" && exit 1

echo ">>> Test finished successfully"
