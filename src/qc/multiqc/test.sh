#!/bin/bash

set -eo pipefail

## VIASH START
meta_temp_dir="/tmp/"
meta_functionality_name="multiqc"
meta_executable="./target/docker/qc/multiqc/multiqc"
meta_resources_dir="resources_test/10x_5k_anticmv/"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT


echo ">>> Running multiqc, outputting to $tmpdir"
# $meta_executable --input "${fastq_files[@]}" --output "$tmpdir"
$meta_executable --input "$meta_resources_dir/fastqc" --output "$tmpdir"

echo ">>> Checking exitcode"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Test failed: non-zero exitcode, was $exit_code" && exit 1
fi

echo ">>> Checking if output was created"
[[ ! -f "$tmpdir/multiqc_report.html" ]] && echo "Output report was not found!" && exit 1

echo ">>> Test finished successfully"
