#!/bin/bash

set -eo pipefail

# TODO: we should turn this into viash components

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=bdrhap_vdj
OUT=resources_test/$ID

# create raw directory
raw_dir="$OUT/raw"
mkdir -p "$raw_dir"

# Check whether seqkit is available
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

sample_ids=(
  RhapVDJDemo-BCR_S1_L001_R1_001 RhapVDJDemo-BCR_S1_L001_R2_001
  RhapVDJDemo-mRNA_S5_L001_R1_001 RhapVDJDemo-mRNA_S5_L001_R2_001
  RhapVDJDemo-TCR_S3_L001_R1_001 RhapVDJDemo-TCR_S3_L001_R2_001
)
fasta_file="$raw_dir/BD_Rhapsody_Immune_Response_Panel_Hs.fasta"

# check whether all raw files have already been created
raw_files_exist=true
for sample_id in "${sample_ids[@]}"; do
  [[ -f "$raw_dir/${sample_id}_subset.fastq.gz" ]] || raw_files_exist=false
done
[[ -f "$fasta_file" ]] || raw_files_exist=false

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/VDJDemo"
if [[ "$raw_files_exist" == false && ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"
    wget "http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-Demo-Data-Inputs/VDJDemo/VDJDemo.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

# subset fastq files
for sample_id in "${sample_ids[@]}"; do
  subset_file="$raw_dir/${sample_id}_subset.fastq.gz"
  if [[ ! -f "$subset_file" ]]; then
    echo "> Processing $sample_id"
    seqkit head -n 300000 "$tar_dir/$sample_id.fastq.gz" | gzip > "$subset_file"
  fi
done

# copy immune panel fasta
if [[ ! -f "$fasta_file" ]]; then
  cp "$tar_dir/BD_Rhapsody_Immune_Response_Panel_Hs.fasta" "$fasta_file"
fi

# run bd rhapsody pipeline
# remove previously generated output, so no outdated files are left behind
rm -rf "$OUT/processed"
processed_dir="$OUT/processed/output_raw"
mkdir -p "$processed_dir"
reads=$(ls "$raw_dir"/RhapVDJDemo-*_S*_L001_R[12]_001_subset.fastq.gz | paste -sd ";" -)
viash run src/mapping/bd_rhapsody/config.vsh.yaml -- \
  --reads "$reads" \
  --targeted_reference "$fasta_file" \
  --vdj_version human \
  --output_dir "$processed_dir"
