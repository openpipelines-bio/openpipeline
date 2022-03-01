#!/usr/bin/env bash

set -ex

test_file="$meta_resources_dir/LICENSE"

echo ">>> Check whether test file exists"
[[ ! -f $test_file ]] && echo "Test file could not be found!" && exit 1

echo ">>> Creating tar.gz..."
tar czvf $test_file.tar.gz $test_file

echo ">>> Check whether tar.gz can be extracted"
./$meta_functionality_name \
   --input "$test_file.tar.gz" \
   --output output_dir

[[ ! -d output_dir ]] && echo "Output directory could not be found!" && exit 1

echo ">>> Test finished successfully"
