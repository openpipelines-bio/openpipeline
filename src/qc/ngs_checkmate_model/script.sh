#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END


if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi

cd NGSCheckMate
python2.7 ncm.py -V -d "$par_input" -bed SNP/SNP_GRCh38_hg38_wChr.bed -O "$par_output"
