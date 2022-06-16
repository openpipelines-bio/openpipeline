#!/bin/bash

bcl-convert \
  --force \
  --bcl-input-directory "$par_input" \
  --output-directory "$par_output" \
  --sample-sheet "$par_sample_sheet" \
  --first-tile-only $par_test_mode
