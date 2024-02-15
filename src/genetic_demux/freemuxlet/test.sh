#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --plp "$meta_resources_dir/demuxafy_test_data/dsc_pileup/dscpileup_out" --output freemuxlet_result/ --out freemux_out --nsample 14

[[ ! -f freemuxlet_result/freemux_out.clust1.samples.gz ]] && echo "Output VCF file could not be found!" && exit 1
[[ ! -f freemuxlet_result/cell_annotation.csv ]] && echo "Output cell type annotation file could not be found!" && exit 1


echo ">>> Test finished successfully"
