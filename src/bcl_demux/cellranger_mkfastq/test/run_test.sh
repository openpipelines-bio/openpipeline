#!/usr/bin/env bash

set -ex

echo ">>> Checking whether output is generated"
./cellranger_mkfastq \
   --input cellranger-tiny-bcl-1.2.0.tar.gz \
   --samplesheet cellranger-tiny-bcl-simple-1.2.0.csv \
   --output tinybcl_output \
   --cores 1 \
   --memory 1

[[ ! -d tinybcl_output ]] && echo "Output file could not be found!" && exit 1

rm -rf ./tinybcl_output 

tar xzf "cellranger-tiny-bcl-1.2.0.tar.gz"

./cellranger_mkfastq \
   --input cellranger-tiny-bcl-1.2.0 \
   --samplesheet cellranger-tiny-bcl-simple-1.2.0.csv \
   --output tinybcl_output_2 \
   --barcodeMismatches 0 \
   --cores 1 \
   --memory 1

[[ ! -d tinybcl_output_2 ]] && echo "Output file could not be found!" && exit 1

echo ">>> Test finished successfully"
