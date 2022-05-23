#!/bin/bash

## VIASH START
par_input='resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/'
par_reference='resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/'
par_output='resources_test/cellranger_tiny_fastq/bam'

## VIASH END

# just to make sure
par_input=`realpath $par_input`
par_reference=`realpath $par_reference`
par_output=`realpath $par_output`

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_resources_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# cd into tempdir
cd "$tmpdir"

# add additional params
extra_params=( )

if [ ! -z "$par_cores" ]; then 
  extra_params+=( "--localcores" "$par_cores" )
fi
if [ ! -z "$par_memory" ]; then 
  extra_params+=( "--localmem" "$par_memory" )
fi
if [ ! -z "$par_expect_cells" ]; then 
  extra_params+=( "--expect-cells" "$par_expect_cells" )
fi
if [ ! -z "$par_chemistry" ]; then 
  extra_params+=( "--chemistry" "$par_chemistry" )
fi

echo "Running cellranger count"


id=myoutput
cellranger count \
  --id "$id" \
  --fastqs "$par_input" \
  --transcriptome "$par_reference" \
  "${extra_params[@]}" \
  --disable-ui \
  --nosecondary

echo "Copying output"
if [ -d "$id/outs/" ]; then
  if [ ! -d "$par_output" ]; then
    mkdir -p "$par_output"
  fi
  mv "$id/outs/"* "$par_output"
fi
