#!/bin/bash

set -eou pipefail

## VIASH START
par_input='resources_test/cellranger_tiny_bcl/bcl'
par_sample_sheet='resources_test/cellranger_tiny_bcl/bcl/sample_sheet.csv'
par_output=foo
par_input=`realpath $par_input`
par_sample_sheet=`realpath $par_sample_sheet`
par_output=`realpath $par_output`
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_name-XXXXXXXX")
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


# add additional params
extra_params=( )

if [ ! -z "$meta_cpus" ]; then 
  extra_params+=( "--localcores=$meta_cpus" )
fi
if [ ! -z "$meta_memory_gb" ]; then 
  # always keep 2gb for the OS itself
  memory_gb=`python -c "print(int('$meta_memory_gb') - 2)"`
  extra_params+=( "--localmem=$memory_gb" )
fi


echo "Running cellranger demux"

id=myoutput

cellranger mkfastq \
  --id "$id" \
  --csv "$par_sample_sheet" \
  --run "$par_input" \
  "${extra_params[@]}" \
  --disable-ui \
  --output-dir "$par_output"

# Move reports to their own output location
if [ ! -z "$par_reports" ]; then
  echo "Moving reports its own location"
  mv "$par_output"/Reports "$par_reports"
else
  echo "Leaving reports alone"
fi
