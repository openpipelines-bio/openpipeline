#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=10x_5k_fixed
OUT="resources_test/$ID"

# create raw directory
raw_dir="$OUT/raw"
mkdir -p "$raw_dir"

# Check whether seqkit is available
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

# check whether reference is available
reference_dir="resources_test/reference_gencodev41_chr1/"
genome_tar="$reference_dir/reference_cellranger.tar.gz"
if [[ ! -f "$genome_tar" ]]; then
    echo "$genome_tar does not exist. Please create the reference genome first"
    exit 1
fi

# dataset page:
# https://www.10xgenomics.com/datasets/Mixture-of-cells-from-mouse-lymph-nodes-and-spleen-stained-with-totalseqc-mouse-universal-cocktail

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"

    # download fastqs and untar
    wget "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Multiplex/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Multiplex_fastqs.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

function seqkit_head {
  input="$1"
  output="$2"
  if [[ ! -f "$output" ]]; then
    echo "> Processing `basename $input`"
    seqkit head -n 200000 "$input" | gzip > "$output"
  fi
}

orig_sample_id="4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex"

seqkit_head "$tar_dir/gex/${orig_sample_id}_ab_S1_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_ab_subset_S1_L001_R1_001.fastq.gz"
seqkit_head "$tar_dir/gex/${orig_sample_id}_ab_S1L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_ab_subset_S1_L001_R2_001.fastq.gz"

seqkit_head "$tar_dir/vdj/${orig_sample_id}_gex_S2_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S2_L001_R1_001.fastq.gz" 
seqkit_head "$tar_dir/vdj/${orig_sample_id}_gex_S2_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S2_L001_R2_001.fastq.gz" 

# download feature reference
feature_ref="$raw_dir/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_feature_reference.csv"
if [[ ! -f "$feature_ref" ]]; then
  wget "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Multiplex/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Multiplex_count_feature_reference.csv" -O "$feature_ref"
fi

# download probe set
probe_set="$raw_dir/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
if [[ ! -f "$probe_set" ]]; then
  wget "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv" -O "$feature_ref"
fi


# Run mapping pipeline
cat > /tmp/params.yaml << HERE
param_list:
- id: "$ID"
  input: "$raw_dir"
  library_id:
    - "${orig_sample_id}_ab_subset"
    - "${orig_sample_id}_gex_subset"
  library_type:
    - "Antigen Capture"
    - "Gene Expression"

probe_set: "$probe_set"
gex_reference: "$genome_tar"
feature_reference: "$feature_ref"
publish_dir: "$OUT/processed"
HERE

nextflow \
  run . \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels_ci.config

# Create h5mu
cat > /tmp/params.yaml << HERE
id: "$ID"
input: "$OUT/processed/$ID.cellranger_multi.output"
publish_dir: "$OUT/"
output: "$orig_sample_id.h5mu"
HERE

nextflow \
  run . \
  -main-script target/nextflow/convert/from_cellranger_multi_to_h5mu/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels_ci.config
