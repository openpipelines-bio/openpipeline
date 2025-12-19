#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=10x_5k_beam
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
# https://www.10xgenomics.com/datasets/5k-human-a0201-b0702-pbmcs-beam-t-2-standard

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/5k_human_A0201_B0702_PBMCs_BEAM_T"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"

    # download fastqs and untar
    wget "https://cf.10xgenomics.com/samples/cell-vdj/7.1.0/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex_fastqs.tar" -O "$tar_dir.tar"
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

orig_sample_id="beamt_human_A0201_B0702_pbmc"

seqkit_head "$tar_dir/gex/${orig_sample_id}_gex_S3_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S3_L001_R1_001.fastq.gz"
seqkit_head "$tar_dir/gex/${orig_sample_id}_gex_S3_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S3_L001_R2_001.fastq.gz"

seqkit_head "$tar_dir/vdj/${orig_sample_id}_vdj_S2_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_vdj_subset_S2_L001_R1_001.fastq.gz" 
seqkit_head "$tar_dir/vdj/${orig_sample_id}_vdj_S2_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_vdj_subset_S2_L001_R2_001.fastq.gz" 

seqkit_head "$tar_dir/antigen_capture/${orig_sample_id}_ag_S1_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_ag_subset_S1_L001_R1_001.fastq.gz" 
seqkit_head "$tar_dir/antigen_capture/${orig_sample_id}_ag_S1_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_ag_subset_S1_L001_R2_001.fastq.gz" 

# download feature reference
feature_ref="$raw_dir/beamt_human_A0201_B0702_pbmc_feature_reference.csv"
if [[ ! -f "$feature_ref" ]]; then
  wget "https://cf.10xgenomics.com/samples/cell-vdj/7.1.0/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex_count_feature_reference.csv" -O "$feature_ref"
fi

# download vdj reference if needed
vdj_ref="$raw_dir/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex_vdj_reference.tar.gz"
if [[ ! -f "$vdj_ref" ]]; then
  wget "https://cf.10xgenomics.com/samples/cell-vdj/7.1.0/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex_vdj_reference.tar.gz" -O "$vdj_ref"
fi

# Run mapping pipeline
# TODO: Also include conversion to h5mu
cat > /tmp/params.yaml << HERE
param_list:
- id: "$ID"
  input: "$raw_dir"
  library_id:
    - "${orig_sample_id}_gex_subset"
    - "${orig_sample_id}_vdj_subset"
    - "${orig_sample_id}_ag_subset"
  library_type:
    - "Gene Expression"
    - "VDJ-T"
    - "Antigen Capture"

gex_reference: "$genome_tar"
feature_reference: "$feature_ref"
vdj_reference: "$vdj_ref"
control_id:
    - negative_control_A0201
    - negative_control_B0702
mhc_allele:
    - "HLA-A*02:01"
    - "HLA-B*07:02"
HERE

nextflow \
  run openpipelines-bio/openpipeline \
  -r 3.0.0 \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$OUT/processed"


nextflow \
  run . \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "${OUT}_v10/processed"

# Create h5mu
cat > /tmp/params.yaml << HERE
id: "$ID"
input: "$OUT/processed/$ID.cellranger_multi.output"
publish_dir: "$OUT/"
output: "*.h5mu"
HERE

nextflow \
  run openpipelines-bio/openpipeline \
  -r 3.0.0 \
  -main-script target/nextflow/convert/from_cellranger_multi_to_h5mu/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels_ci.config

mv "$OUT/0.h5mu" "$OUT/${orig_sample_id}.h5mu"