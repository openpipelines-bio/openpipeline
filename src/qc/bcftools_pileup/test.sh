#!/bin/bash

set -eo pipefail

## VIASH START
meta_resources_dir="./resources_test"
meta_executable="./target/docker/qc/bcftools_pileup/bcftools_pileup"
## VIASH END

# Add paths. How to add reference? Also should we do a check for the reference?
"$meta_executable" \
  --reference "" \
  --input "" \
  --output "" 

if [ ! -f "./star_reference_test/Genome" ]; then
    echo "Genome file could not be found in the output directory";
    exit 1
fi