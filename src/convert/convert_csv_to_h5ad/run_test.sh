#!/usr/bin/env bash

set -ex


./convert_csv_to_h5ad --input CS0000007_subsample_LI00080.csv.gz --output output.gz.h5ad

[[ ! -f output.gz.h5ad ]] && echo "Output file could not be found!" && exit 1

gunzip CS0000007_subsample_LI00080.csv

./convert_csv_to_h5ad --input CS0000007_subsample_LI00080.csv --output output.h5ad

[[ ! -f output.h5ad ]] && echo "Output file could not be found!" && exit 1

echo ">>> Test finished successfully"

