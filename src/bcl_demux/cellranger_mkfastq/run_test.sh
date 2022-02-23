#!/usr/bin/env bash

set -ex

data_dir="$resources_dir/cellranger_tiny_bcl_1.2.0"

echo ">>> Check whether component with BCL archives"
./$meta_functionality_name \
   --input "$data_dir/cellranger-tiny-bcl-1.2.0.tar.gz" \
   --samplesheet "$data_dir/cellranger-tiny-bcl-simple-1.2.0.csv" \
   --output tinybcl_output_1 \
   --cores 1 \
   --memory 1

[[ ! -d tinybcl_output_1 ]] && echo "Output file could not be found!" && exit 1



echo ">>> Check whether component with BCL directories"
tar xzf "$data_dir/cellranger-tiny-bcl-1.2.0.tar.gz"

./$meta_functionality_name \
   --input ./cellranger-tiny-bcl-1.2.0 \
   --samplesheet "$data_dir/cellranger-tiny-bcl-simple-1.2.0.csv" \
   --output tinybcl_output_2 \
   --barcodeMismatches 0 \
   --cores 1 \
   --memory 1
[[ ! -d tinybcl_output_2 ]] && echo "Output file could not be found!" && exit 1



echo ">>> Test finished successfully"
