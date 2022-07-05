#!/bin/bash

# settings
ID=bd_rhapsody_wta_5kjrt
OUT=resources_test/$ID
raw_dir="$OUT/raw"
extract_dir="$OUT/12WTA-ABC-SMK-EB-5kJRT_raw" # to be removed before uploading

# create output directory
mkdir -p "$raw_dir"
mkdir -p "$extract_dir"

# Check whether seqkit is available
# TODO: we should turn this into viash components
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

# download 
n_cores=30
wget "http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-Demo-Data-Inputs/12WTA-ABC-SMK-EB-5kJRT.tar" -O "$raw_dir/12WTA-ABC-SMK-EB-5kJRT.tar"
tar -xvf  "$raw_dir/12WTA-ABC-SMK-EB-5kJRT.tar" -C "$OUT/12WTA-ABC-SMK-EB-5kJRT_raw" --strip-components=1
seqkit sample -s100 -p 0.01 -j $n_cores "$extract_dir/12SMK_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12SMK_S1_L432_R1_001.fastq.gz"
seqkit sample -s100 -p 0.01 -j $n_cores "$extract_dir/12SMK_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12SMK_S1_L432_R2_001.fastq.gz"
seqkit sample -s100 -p 0.01 -j $n_cores "$extract_dir/12ABC_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12ABC_S1_L432_R1_001.fastq.gz"
seqkit sample -s100 -p 0.01 -j $n_cores "$extract_dir/12ABC_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12ABC_S1_L432_R2_001.fastq.gz"
seqkit sample -s100 -p 0.01 -j $n_cores "$extract_dir/12WTA_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12WTA_S1_L432_R1_001.fastq.gz"
seqkit sample -s100 -p 0.01 -j $n_cores "$extract_dir/12WTA_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12WTA_S1_L432_R2_001.fastq.gz"

cp "$raw_dir/12WTA-ABC-SMK-EB-5kJRT/BDAbSeq_ImmuneDiscoveryPanel.fasta" "$raw_dir
rm -r "$raw_dir/12WTA-ABC-SMK-EB-5kJRT"
rm "$raw_dir/12WTA-ABC-SMK-EB-5kJRT.tar"

# # process raw files
# processed_dir="$OUT/processed"
# mkdir -p "$processed_dir"

# # run bdrhap
# bdrhap_out="$processed_dir/bdrhap_out/"
# mkdir -p "$bdrhap_out"

# target/docker/mapping/bd_rhapsody_wta/bd_rhapsody_wta \
#   --input "$raw_dir/12SMK_S1_L432_R1_001.fastq.gz" \
#   --input "$raw_dir/12SMK_S1_L432_R2_001.fastq.gz" \
#   --reference_genome "$raw_dir/GRCh38_primary_assembly_genome_chr1.tar.gz" \
#   --transcriptome_annotation "$raw_dir/gencode_v40_annotation_chr1.gtf" \
#   --output "$bdrhap_out"



