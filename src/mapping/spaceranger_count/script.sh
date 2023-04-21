#!/bin/bash

set -eo pipefail

## VIASH START
par_input='resources_test/spaceranger_tiny_fastq/spaceranger_tiny_fastq/'
par_reference='resources_test/spaceranger_tiny_fastq/spaceranger_tiny_ref/'
par_image='sample345.tiff'
par_slide='V19J01-123'
par_area='A1'
par_unknown_slide=false
par_output='resources_test/spaceranger_tiny_fastq/bam'
## VIASH END

# set image type
image='image'

# just to make sure paths are absolute
par_reference=$(realpath "$par_reference")
par_output=$(realpath "$par_output")

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# process inputs
# for every fastq file found, make a symlink into the tempdir
fastq_dir="$tmpdir/fastqs"
mkdir -p "$fastq_dir"
IFS=";"
for var in $par_input; do
  unset IFS
  abs_path=$(realpath "$var")
  if [[ -d "$abs_path" ]]; then
    find "$abs_path" -name *.fastq.gz -exec ln -s {} "$fastq_dir" \;
  else
    ln -s "$abs_path" "$fastq_dir"
  fi
done

# process reference
if file "$par_reference" | grep -q 'gzip compressed data'; then
  echo "Untarring genome"
  reference_dir="$tmpdir/fastqs"
  mkdir -p "$reference_dir"
  tar -xvf "$par_reference" -C "$reference_dir" --strip-components=1
  par_reference="$reference_dir"
fi

# cd into tempdir
cd "$tmpdir"

# add additional params
extra_params=( )

if [[ ! -z "$meta_cpus" ]]; then 
  extra_params+=( "--localcores=$meta_cpus" )
fi
if [[ ! -z "$meta_memory_gb" ]]; then 
  # always keep 2gb for the OS itself
  memory_gb=$(python -c "print(int('$meta_memory_gb') - 2)")
  extra_params+=( "--localmem=$memory_gb" )
fi
if [[ ! -z "$par_slide" ]]; then 
  extra_params+=( "--slide=$par_slide" )
fi
if [[ ! -z "$par_area" ]]; then 
  extra_params+=( "--area=$par_area" )
fi
if [[ ! -z "$par_slidefile" ]]; then
  extra_params+=( "--slidefile" )
fi
if [[ ! -z "$par_unknown_slide" ]]; then
  extra_params+=( "--unknown_slide" )
fi

# Set the variable "image" based on the value of "par_image_type"
if [[ "${par_image_type}" == "darkimage" ]]; then
    image="darkimage"
elif [[ "${par_image_type}" == "colorizedimage" ]]; then
    image="colorizedimage"
else
    echo "Invalid value for par_image_type. Exiting."
    exit 1
fi

echo "Running spaceranger count"

id=myoutput
spaceranger count \
  --id "$id" \
  --fastqs "$fastq_dir" \
  --transcriptome "$par_reference" \
  --"$image" "$par_image" \
  "${extra_params[@]}" \

echo "Copying output"
if [[ -d "$id/outs/" ]]; then
  if [[ ! -d "$par_output" ]]; then
    mkdir -p "$par_output"
  fi
  mv "$id/outs/"* "$par_output"
fi
