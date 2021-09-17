#!/bin/bash

if [ $par_memory == "auto" ]; then
  par_memory=0
fi
if [ $par_cores == "auto" ]; then
  par_cores=0
fi

# making all IO absolute
par_input=$(realpath --no-symlinks "$par_input")
par_samplesheet=$(realpath --no-symlinks "$par_samplesheet")
par_output=$(realpath --no-symlinks "$par_output")

# if par_input is a folder, untar first
if [ ! -d "$par_input" ]; then
  tmpdir=$(mktemp -d)
  function clean_up {
    rm -rf "$tmpdir"
  }
  trap clean_up EXIT
  
  echo "Assuming input is a tar.gz, untarring"
  tar xzf "$par_input" -C "$tmpdir" --strip-components=1
  input_dir="$tmpdir"
else
  input_dir="$par_input"
fi

echo "Running cellranger"
tmpwddir=$(mktemp -d)
function clean_up {
  rm -rf "$tmpwddir"
}
trap clean_up EXIT

cd "$tmpwddir"
cellranger mkfastq \
  --run "$input_dir" \
  --id myoutput \
  --csv "$par_samplesheet" \
  --disable-ui \
  --localmem=$par_memory \
  --barcode-mismatches=$par_barcode_mismatches

echo "Moving output to output directory"
mv "myoutput/outs/fastq_path" "$par_output"

# todo: could interpret Stats/Stats.json to see which output dirs were created
