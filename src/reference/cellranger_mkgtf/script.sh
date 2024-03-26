#!/bin/bash

set -eo pipefail

## VIASH START
par_input_gtf="resources_test/reference_gencodev41_chr1/reference.gtf.gz"
par_output_gtf="gencode_v41_filtered.gtf.gz"
par_attribute="gene_type:transcribed_unprocessed_pseudogene"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# just to make sure
par_input_gtf=`realpath $par_input_gtf`
par_output_gtf=`realpath $par_output_gtf`

echo "> Unzipping input files"
unpigz -c "$par_input_gtf" > "$tmpdir/input_gtf.gtf"

echo "> Building gtf"
cd "$tmpdir"
cellranger mkgtf \
  "$tmpdir/input_gtf.gtf" \
  "$tmpdir/output.gtf" \
  --attribute="${par_attribute}"

echo "> Creating archive"
pigz -k "$tmpdir/output.gtf"
mv "$tmpdir/output.gtf.gz" "$par_output_gtf"