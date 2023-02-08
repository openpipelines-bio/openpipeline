#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

# just to make sure paths are absolute
par_reference=`realpath $par_reference`

samtools mpileup -I -uf "$par_reference" -l NGSCheckMate/SNP/SNP_GRCh38_hg38_wChr.bed "$par_input" | bcftools view -cg - > "$par_output"
