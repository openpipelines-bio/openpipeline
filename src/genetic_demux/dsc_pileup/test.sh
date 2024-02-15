#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --sam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --vcf "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf" \
    --output dscpileup_result/ \
    --out out

[[ ! -f dscpileup_result/out.plp.gz ]] && echo "Output result file could not be found!" && exit 1

echo ">>> Test finished successfully"
