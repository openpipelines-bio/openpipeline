#!/bin/bash

## VIASH START
par_input="input.bam"
par_transcriptome="reference.gtf"
par_barcode="barcodes.tsv"
par_output="output.loom"
## VIASH END

extra_params=( )

if [ ! -z "$par_barcode" ]; then 
  extra_params+=( "--bcfile=$par_barcode" )
fi

velocyto run \
  "$par_input" \
  "$par_transcriptome" \
  "${extra_params[@]}" \
  --outputfolder `dirname $par_output` \
  --sampleid `basename $par_output .loom` \
  --samtools-threads 30 \
  --samtools-memory 3500


