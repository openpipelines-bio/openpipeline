#!/bin/bash

set -eo pipefail

# Handle reports stored separate
reports_line=""
$par_separate_reports && reports_line="--reports-dir $par_reports"

# Handle the boolean flag
ignore=""
$par_ignore_missing && ignore="--ignore-missing-control --ignore-missing-bcl --ignore-missing-filter"

# Run the actual command
bcl2fastq \
  --runfolder-dir "$par_input" \
  --sample-sheet "$par_sample_sheet" \
  --output-dir "$par_output" \
  $reports_line \
  $ignore
