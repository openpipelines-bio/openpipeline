#!/bin/bash

set -eo pipefail

input_dir="$meta_resources_dir/cellranger_count_output"

assert_identical() {
  local expected="$1"
  local actual="$2"
  if [ ! -f "$actual" ]; then
    echo "Expected output file '$actual' does not exist." && exit 1
  fi
  if ! cmp -s "$expected" "$actual"; then
    echo "Output file '$actual' differs from input '$expected'." && exit 1
  fi
}

echo ">>> Test 1: split all output files"
mkdir -p test_all
$meta_executable \
  --input "$input_dir" \
  --filtered_h5 test_all/filtered.h5 \
  --metrics_summary test_all/metrics.csv \
  --molecule_info test_all/molecule_info.h5 \
  --bam test_all/out.bam \
  --bai test_all/out.bam.bai \
  --raw_h5 test_all/raw.h5

assert_identical "$input_dir/filtered_feature_bc_matrix.h5" test_all/filtered.h5
assert_identical "$input_dir/metrics_summary.csv" test_all/metrics.csv
assert_identical "$input_dir/molecule_info.h5" test_all/molecule_info.h5
assert_identical "$input_dir/possorted_genome_bam.bam" test_all/out.bam
assert_identical "$input_dir/possorted_genome_bam.bam.bai" test_all/out.bam.bai
assert_identical "$input_dir/raw_feature_bc_matrix.h5" test_all/raw.h5

echo ">>> Test 2: only requested outputs are created"
mkdir -p test_subset
$meta_executable \
  --input "$input_dir" \
  --filtered_h5 test_subset/filtered.h5 \
  --metrics_summary test_subset/metrics.csv

assert_identical "$input_dir/filtered_feature_bc_matrix.h5" test_subset/filtered.h5
assert_identical "$input_dir/metrics_summary.csv" test_subset/metrics.csv
n_files=$(find test_subset -type f | wc -l)
if [ "$n_files" -ne 2 ]; then
  echo "Expected 2 output files, found $n_files:" && ls -la test_subset && exit 1
fi

echo ">>> Test 3: input directory containing only a subset of the files"
# e.g. Cell Ranger run with --create-bam false: no BAM/BAI in the output directory
mkdir -p partial_input test_partial
cp "$input_dir/filtered_feature_bc_matrix.h5" partial_input/
cp "$input_dir/raw_feature_bc_matrix.h5" partial_input/
cp "$input_dir/metrics_summary.csv" partial_input/
$meta_executable \
  --input partial_input \
  --filtered_h5 test_partial/filtered.h5 \
  --raw_h5 test_partial/raw.h5 \
  --metrics_summary test_partial/metrics.csv

assert_identical "$input_dir/filtered_feature_bc_matrix.h5" test_partial/filtered.h5
assert_identical "$input_dir/raw_feature_bc_matrix.h5" test_partial/raw.h5
assert_identical "$input_dir/metrics_summary.csv" test_partial/metrics.csv
n_files=$(find test_partial -type f | wc -l)
if [ "$n_files" -ne 3 ]; then
  echo "Expected 3 output files, found $n_files:" && ls -la test_partial && exit 1
fi

echo ">>> All tests finished successfully"
