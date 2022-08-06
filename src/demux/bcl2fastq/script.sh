#!/bin/bash

ignore=""

$par_ignore_missing && ignore="--ignore-missing-control --ignore-missing-bcl --ignore-missing-filter"

bcl2fastq \
  --runfolder-dir "$par_input" \
  --sample-sheet "$par_sample_sheet" \
  --output-dir "$par_output" \
  $ignore
