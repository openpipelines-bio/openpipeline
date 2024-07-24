#!/usr/bin/env bash

set -ex

INPUT_FILE="LICENSE"
OUTPUT_DIR="output/files/"
OUTPUT_FILE="${OUTPUT_DIR}${INPUT_FILE}"

echo ">>> Check whether test file exists"
[[ ! -f ${INPUT_FILE} ]] && echo "Test file could not be found!" && exit 1

echo ">>> Creating tar.gz..."
tar czvf ${INPUT_FILE}.tar.gz ${INPUT_FILE}

echo ">>> Check whether tar.gz can be extracted"
./$meta_executable \
   --input "${INPUT_FILE}.tar.gz" \
   --output "$OUTPUT_DIR"

echo ">>> Check whether extracted file exists"
[[ ! -f $OUTPUT_FILE ]] && echo "Output file could not be found!" && exit 1

echo ">>> Check whether input and output file are the same"
cmp $INPUT_FILE $OUTPUT_FILE || (echo "Input and output files are different!" && exit 1)

echo ">>> Test finished successfully"
