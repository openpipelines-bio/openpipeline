#!/bin/bash

# settings
ID=bd_rhapsody_wta_1_10_test
OUT=resources_test/$ID
DIR="$OUT"
S3DIR=$(echo "$DIR" | sed 's#resources_test#s3://openpipelines-data#')
raw_dir="$OUT/raw"


mkdir -p "$raw_dir"
wget http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-Demo-Data-Inputs/12WTA-ABC-SMK-EB-5kJRT.tar -O "$raw_dir/12WTA-ABC-SMK-EB-5kJRT.tar"
mkdir "$raw_dir/12WTA-ABC-SMK-EB-5kJRT"
tar -xvf  "$raw_dir/12WTA-ABC-SMK-EB-5kJRT.tar" -C "$raw_dir/12WTA-ABC-SMK-EB-5kJRT" --strip-components=1
seqkit sample -s100 -p 0.01 "$raw_dir/12WTA-ABC-SMK-EB-5kJRT/12SMK_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12SMK_S1_L432_R1_001.fastq.gz"
seqkit sample -s100 -p 0.01 "$raw_dir/12WTA-ABC-SMK-EB-5kJRT/12SMK_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12SMK_S1_L432_R2_001.fastq.gz"
seqkit sample -s100 -p 0.01 "$raw_dir/12WTA-ABC-SMK-EB-5kJRT/12ABC_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12ABC_S1_L432_R1_001.fastq.gz"
seqkit sample -s100 -p 0.01 "$raw_dir/12WTA-ABC-SMK-EB-5kJRT/12ABC_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12ABC_S1_L432_R2_001.fastq.gz"
seqkit sample -s100 -p 0.01 "$raw_dir/12WTA-ABC-SMK-EB-5kJRT/12WTA_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12WTA_S1_L432_R1_001.fastq.gz"
seqkit sample -s100 -p 0.01 "$raw_dir/12WTA-ABC-SMK-EB-5kJRT/12WTA_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12WTA_S1_L432_R1_001.fastq.gz"

cp "$raw_dir/12WTA-ABC-SMK-EB-5kJRT/BDAbSeq_ImmuneDiscoveryPanel.fasta" "$raw_dir
rm -r "$raw_dir/12WTA-ABC-SMK-EB-5kJRT"
rm "$raw_dir/12WTA-ABC-SMK-EB-5kJRT.tar"

# # process raw files
# processed_dir="$OUT/processed"
# mkdir -p "$processed_dir"

# # run bdrhap
# bdrhap_out="$processed_dir/bdrhap_out/"
# mkdir -p "$bdrhap_out"

# target/docker/mapping/bd_rhapsody_wta_1_10_1/bd_rhapsody_wta_1_10_1 \
#   --input "$raw_dir/12SMK_S1_L432_R1_001.fastq.gz" \
#   --input "$raw_dir/12SMK_S1_L432_R2_001.fastq.gz" \
#   --reference_genome "$raw_dir/GRCh38_primary_assembly_genome_chr1.tar.gz" \
#   --transcriptome_annotation "$raw_dir/gencode_v40_annotation_chr1.gtf" \
#   --output "$bdrhap_out"




# aws s3 sync --profile xxx "$DIR" "$S3DIR"
