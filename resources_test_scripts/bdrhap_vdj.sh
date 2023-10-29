#!/bin/bash

set -eo pipefail

# TODO: we should turn this into viash components

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=bdrhap_vdj
OUT=resources_test/$ID
n_threads=30

# create raw directory
raw_dir="$OUT/raw"
mkdir -p "$raw_dir"

# Check whether seqkit is available
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/VDJDemo"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"
    wget "http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-Demo-Data-Inputs/VDJDemo/VDJDemo.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

# subset fastq files
for sample_id in RhapVDJDemo-BCR_S1_L001_R1_001 RhapVDJDemo-BCR_S1_L001_R2_001 RhapVDJDemo-mRNA_S5_L001_R1_001 RhapVDJDemo-mRNA_S5_L001_R2_001 RhapVDJDemo-TCR_S3_L001_R1_001 RhapVDJDemo-TCR_S3_L001_R2_001; do
  subset_file="$raw_dir/${sample_id}_subset.fastq.gz"
  if [[ ! -f "$subset_file" ]]; then
  echo "> Processing $sample_id"
    seqkit head -n 300000 "$tar_dir/$sample_id.fastq.gz" | gzip > "$subset_file"
  fi
  unset subset_file
done

# copy immune panel fasta
fasta_file="$raw_dir/BD_Rhapsody_Immune_Response_Panel_Hs.fasta"
if [[ ! -f "$fasta_file" ]]; then
  cp "$tar_dir/BD_Rhapsody_Immune_Response_Panel_Hs.fasta" "$fasta_file"
fi

# create params file
cat > /tmp/params.yaml << HERE
param_list:
- id: "targeted_vdj"
  input: "$raw_dir/RhapVDJDemo-*_S*_L001_R[12]_001_subset.fastq.gz"
mode: targeted
reference: "$fasta_file"
publish_dir: "$OUT/processed"
putative_cell_call: "mRNA"
vdj_version: human
HERE

# run bd rhapsody pipeline
nextflow \
  run . \
  -main-script src/workflows/ingestion/bd_rhapsody/main.nf \
  -resume \
  -profile docker,mount_temp \
  -with-trace work/trace.txt \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels.config \
  -c src/workflows/utils/errorstrat_ignore.config