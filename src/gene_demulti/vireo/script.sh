#!/bin/bash

set -eo pipefail

# Unset flags if they equal 'false'
[[ "$par_noDoublet" == "false" ]] && unset par_noDoublet
[[ "$par_forceLearnGT" == "false" ]] && unset par_forceLearnGT
[[ "$par_ASEmode" == "false" ]] && unset par_ASEmode
[[ "$par_noPlot" == "false" ]] && unset par_noPlot
[[ "$par_callAmbientRNAs" == "false" ]] && unset par_callAmbientRNAs

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

vireo \
  --cellData $par_cellData \
  --nDonor $par_nDonor \
  --genoTag $par_genoTag \
  --nInit $par_nInit \
  --extraDonor $par_extraDonor \
  --nproc $par_nproc \
  --out "${par_output}" \
  ${par_vartrixData:+--vatrixData $par_vartrixData} \
  ${par_donorFile:+--donorFile $par_donorFile} \
  ${par_noDoublet:+--noDoublet} \
  ${par_extraDonorMode:+--extraDonorMode $par_extraDonorMode} \
  ${par_forceLearnGT:+--forceLearnGT} \
  ${par_ASEmode:+--ASEmode} \
  ${par_noPlot:+--noPlot} \
  ${par_randSeed:+--randSeed $par_randSeed} \
  ${par_cellRange:+--cellRange $par_cellRange} \
  ${par_callAmbientRNAs:+--callAmbientRNAs}

cut -d$'\t' -f 1-2 "${par_output}/donor_ids.tsv" > "${par_output}/assignment.tsv"
