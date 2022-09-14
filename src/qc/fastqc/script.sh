#!/bin/bash

set -eo pipefail

mkdir -p "$par_output"

if [ "$par_mode" == "dir" ]; then
  par_input="$par_input/*.fastq.gz"
fi

eval fastqc -o "$par_output" "$par_input"
