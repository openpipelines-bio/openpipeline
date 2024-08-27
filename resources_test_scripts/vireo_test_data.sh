#!/bin/bash

set -eo pipefail


# settings
ID=vireo_test_data
OUT=resources_test/$ID
DIR="$OUT"

mkdir -p "$OUT"
cd "$OUT"
# download vireo tutorial dataset
wget https://github.com/single-cell-genetics/vireo/raw/master/data/cells.cellSNP.vcf.gz
