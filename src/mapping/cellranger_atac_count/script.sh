#!/bin/bash

set -eo pipefail

## VIASH START
par_input='resources_test/cellranger_atac_tiny_bcl/fastqs/'
par_reference='resources_test/reference_gencodev41_chr1/'
par_output='resources_test/cellranger_atac_tiny_bcl/bam'
## VIASH END

# just to make sure paths are absolute
par_reference=`realpath $par_reference`
par_output=`realpath $par_output`

echo "Creating temporary directory"
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# process inputs
# for every fastq file found, make a symlink into the tempdir
echo "Locating fastqs"
fastq_dir="$tmpdir/fastqs"
mkdir -p "$fastq_dir"
IFS=";"
for var in $par_input; do
  unset IFS
  abs_path=`realpath $var`
  if [ -d "$abs_path" ]; then
    find "$abs_path" -name *.fastq.gz -exec ln -s {} "$fastq_dir" \;
  else
    ln -s "$abs_path" "$fastq_dir"
  fi
done

echo "fastq_dir content: $(ls $fastq_dir)"

echo "Processing reference"
# process reference
if file $par_reference | grep -q 'gzip compressed data'; then
  echo "Untarring genome"
  reference_dir="$tmpdir/fastqs"
  mkdir -p "$reference_dir"
  tar -xvf "$par_reference" -C "$reference_dir" --strip-components=1
  par_reference="$reference_dir"
fi

# cd into tempdir
cd "$tmpdir"

if [ ! -z "$meta_memory_gb" ]; then 
  # always keep 2gb for the OS itself
  memory_gb=`python -c "print(int('$meta_memory_gb') - 2)"`
fi

echo "Running cellranger-atac count"

id=myoutput
cellranger-atac count \
  --id "$id" \
  --fastqs "$fastq_dir" \
  --reference "$par_reference" \
  --dim-reduce "$par_dim_reduce" \
  --description "$par_description" \
  ${par_lanes:+--lanes=${par_lanes[*]}} \
  ${par_force_cells:+--force-cells=$par_force_cells} \
  ${par_subsample_rate:+--subsample-rate=$par_subsample_rate} \
  ${memory_gb:+--localmem=$memory_gb} \
  ${meta_cpus:+--localcores=$meta_cpus} \
  ${par_lanes:+--lanes=${par_lanes[*]}}

echo "Copying output"
if [ -d "$id/outs/" ]; then
  if [ ! -d "$par_output" ]; then
    mkdir -p "$par_output"
  fi
  mv "$id/outs/"* "$par_output"
fi
