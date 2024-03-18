#!/bin/bash

set -eo pipefail

## VIASH START
par_input='resources_test/cellranger_atac_tiny_bcl/bcl'
par_csv='resources_test/cellranger_atac_tiny_bcl/bcl/sample_sheet.csv'
par_output=foo

par_input=`realpath $par_input`
par_sample_sheet=`realpath $par_csv`
par_output=`realpath $par_output`
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# if par_input not is a folder, untar first
if [ ! -d "$par_input" ]; then  
  echo "Assuming input is a tar.gz, untarring"
  input_dir="$tmpdir/bcl"
  mkdir -p "$input_dir"
  tar -xzf "$par_input" -C "$input_dir" --strip-components=1
else
  input_dir="$par_input"
fi


if [ ! -z "$meta_memory_gb" ]; then 
  # always keep 2gb for the OS itself
  memory_gb=`python -c "print(int('$meta_memory_gb') - 2)"`
fi


echo "Running cellranger-atac mkfastq"

id=myoutput

IFS=","
cellranger-atac mkfastq \
  --id "$id" \
  --csv "$par_csv" \
  --run "$par_input" \
  --disable-ui \
  --output-dir "$par_output" \
  ${meta_cpus:+--localcores=$meta_cpus} \
  ${memory_gb:+--localmem=$memory_gb} \
  ${par_lanes:+--lanes=${par_lanes[*]}} \
  ${par_use_bases_mask:+--use-bases-mask=${par_use_bases_mask[*]} \
  ${par_delete_undetermined:+--delete-undetermined} \
  ${par_barcode_mismatches:+--barcode-mismatches=$par_barcode_mismatches}
unset IFS

# Move reports to their own output location
if [ ! -z "$par_reports" ]; then
  echo "Moving reports its own location"
  mv "$par_output"/Reports "$par_reports"
else
  echo "Leaving reports alone"
fi
