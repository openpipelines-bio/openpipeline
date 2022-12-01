#!/bin/bash

set -exo pipefail

extra_params=()

# Handle reports stored separate
if [ ! -z "$par_reports" ]; then
  extra_params+=("--reports-dir" "$par_reports")
fi

# Handle the boolean flag
if [ "$par_ignore_missing" == "true" ]; then
  extra_params+=("--ignore-missing-control" "--ignore-missing-bcl" "--ignore-missing-filter")
fi

# Run the actual command
bcl2fastq \
  --runfolder-dir "$par_input" \
  --sample-sheet "$par_sample_sheet" \
  --output-dir "$par_output" \
  "${extra_params[@]}"
