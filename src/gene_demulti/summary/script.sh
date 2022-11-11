#!/bin/bash

set -eo pipefail
if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi
Rscript summary/compareParam.R --file $par_file --barcode $par_barcode --outdir $par_output