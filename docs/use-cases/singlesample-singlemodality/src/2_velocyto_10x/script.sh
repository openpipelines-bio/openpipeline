#!/bin/bash

## VIASH START
par_input="fastq_dir"
par_transcriptome="reference_dir"
par_output="output"
## VIASH END


# defaults based on :
# https://github.com/velocyto-team/velocyto.py/blob/112a21601909f371ca061cd86d2592874a0e399e/velocyto/commands/run10x.py#L84

bam_file="$par_input/possorted_genome_bam.bam"

if [ ! -f "$bam_file" ]; then
  echo "Could not locate BAM file at: $bam_file"
  exit 1
fi

barcode_file="$par_input/filtered_feature_bc_matrix/barcodes.tsv.gz"
if [ ! -f "$barcode_file" ]; then
  echo "Could not locate barcodes.tsv.gz file at: $barcode_file"
  exit 1
fi

gtf_file="$par_transcriptome/genes/genes.gtf"
if [ ! -f "$gtf_file" ]; then
  echo "Could not locate GTF file at: $gtf_file"
  exit 1
fi

# autodetect number of threads
num_threads=`nproc --all`

# autodetect 85% of available memory
prec=0
read -r _ freemem _ <<< "$(grep --fixed-strings 'MemFree' /proc/meminfo)"
amnt_mem=`bc <<< "scale=${prec:-3};${freemem}/1024/$num_threads/0.85"`

velocyto run \
  "$bam_file" \
  "$gtf_file" \
  --bcfile "$barcode_file" \
  --outputfolder `dirname $par_output` \
  --sampleid `basename $par_output .loom` \
  --samtools-threads $num_threads \
  --samtools-memory $amnt_mem


