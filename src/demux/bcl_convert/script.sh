#!/bin/bash

bcl-convert \
  --force \
  --bcl-input-directory "$par_input" \
  --output-directory "$par_output" \
  --sample-sheet "$par_sample_sheet" \
  --first-tile-only $par_test_mode

if [ "$par_separate_reports" == true ]; then
  echo "Moving reports to its own location"
  mv "$par_output"/Reports "$par_reports"
else
  echo "Leaving reports alone"
fi
