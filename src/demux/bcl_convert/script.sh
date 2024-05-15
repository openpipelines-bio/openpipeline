#!/bin/bash

set -eo pipefail

## VIASH START
par_input=work/fc/34b01dbb67178188ce8571b7c5459e/bcl2
par_output=work/fc/34b01dbb67178188ce8571b7c5459e/foo
par_sample_sheet=work/fc/34b01dbb67178188ce8571b7c5459e/sample_sheet.csv
par_test_mode=false
## VIASH END

[ -d "$par_output" ] || mkdir -p "$par_output"

bcl-convert \
  --force \
  --bcl-input-directory "$par_input" \
  --output-directory "$par_output" \
  --sample-sheet "$par_sample_sheet" \
  --first-tile-only "$par_test_mode" \
  --strict-mode "$par_strict_mode" \
  ${par_no_lane_splitting:+--no-lane-splitting "$par_no_lane_splitting"} \
  ${par_tiles:+--tiles $par_tiles} \
  ${par_exclude_tiles:+--exclude-tiles $par_exclude_tiles} 
  

if [ ! -z "$par_reports" ]; then
  echo "Moving reports to its own location"
  mv "$par_output"/Reports "$par_reports"
else
  echo "Leaving reports alone"
fi