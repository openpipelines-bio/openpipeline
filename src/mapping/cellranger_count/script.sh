#!/bin/bash

set -eo pipefail

## VIASH START
par_input='resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/'
par_reference='resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/'
par_output='resources_test/cellranger_tiny_fastq/bam'
par_chemistry="auto"
par_expect_cells="3000"
par_secondary_analysis="false"
## VIASH END

# just to make sure paths are absolute
par_reference=`realpath $par_reference`
par_output=`realpath $par_output`

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
  abs_path=`realpath $var`
  if [ -d "$abs_path" ]; then
    find "$abs_path" -name *.fastq.gz -exec ln -s {} "$fastq_dir" \;
  else
    ln -s "$abs_path" "$fastq_dir"
  fi
done

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
if [ ! -z "$par_expect_cells" ]; then 
  extra_params+=( "--expect-cells=$par_expect_cells" )
fi
if [ ! -z "$par_chemistry" ]; then 
  extra_params+=( "--chemistry=$par_chemistry" )
fi
if [ "$par_secondary_analysis" == "false" ]; then
  extra_params+=( "--nosecondary" )
fi
if [ "$par_generate_bam" == "false" ]; then
  extra_params+=( "--no-bam" )
fi
echo "Running cellranger count"


id=myoutput
cellranger count \
  --id "$id" \
  --fastqs "$fastq_dir" \
  --transcriptome "$par_reference" \
  --include-introns "$par_include_introns" \
  "${extra_params[@]}" \
  --disable-ui \

echo "Copying output"
if [ -d "$id/outs/" ]; then
  if [ ! -d "$par_output" ]; then
    mkdir -p "$par_output"
  fi
  mv "$id/outs/"* "$par_output"
fi
