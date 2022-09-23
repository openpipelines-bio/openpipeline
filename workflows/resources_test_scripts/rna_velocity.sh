#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=rna_velocity
OUT=resources_test/$ID


# create raw directory
velocyto_dir="$OUT/velocyto"
mkdir -p "$velocyto_dir"

########################################################
# Create a compatible BAM file from BD Rhapsody Output #
########################################################

bd_rhap_wta_bam="resources_test/bdrhap_5kjrt/processed/SMK.bd_rhapsody.output/Combined_sample_final.BAM"

if [[ ! -f "$bd_rhap_wta_bam" ]]; then
    echo "$bd_rhap_wta_bam does not exist. Please generate BD Rhapsody test data first."
    exit 1
fi

echo "> Converting BD Rhapsody barcode tags."
./bin/viash run ./src/convert/from_10x_to_bd_molecular_barcode_tags/config.vsh.yaml -- \
  -i "$bd_rhap_wta_bam" \
  -o "$velocyto_dir/compatible_bd_input.bam" \
  --bam \
  -t 4

echo "> Creating barcodes file."
samtools view -@4 "$velocyto_dir/compatible_bd_input.bam" | \
  grep -oP "(?<=CB:Z:)\S+" | sort | uniq | head > "$velocyto_dir/barcodes.txt"

###########################################################
# Process Tiny Fast Fastq dataset from 10X to create      #
# input data for convert/from_velocyto_to_h5mu compontent #
###########################################################

mkdir "$OUT/velocyto_processed"

gtf="resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz"
gtf_decompressed="$OUT/velocyto_processed/genes.gtf"
bam="resources_test/cellranger_tiny_fastq/bam/possorted_genome_bam.bam"

echo "> Decompressing gene annotation"
gzip -k -d -c "$gtf" > "$gtf_decompressed"

echo "> Processing 10x dataset"
./bin/viash run ./src/projection/velocyto/config.vsh.yaml -- \
  -i "$bam" \
  -o "$OUT/velocyto_processed/cellranger_tiny.loom" \
  --transcriptome "$gtf_decompressed"

rm "$gtf_decompressed"