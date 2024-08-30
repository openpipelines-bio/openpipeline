#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=10x_5k_lung_crispr
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
# https://www.10xgenomics.com/resources/datasets/5-k-a-549-lung-carcinoma-cells-no-treatment-transduced-with-a-crispr-pool-3-1-standard-6-0-0

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"

    # download fastqs and untar
    wget "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_fastqs.tar" -O "$tar_dir.tar"
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

orig_sample_id="SC3_v3_NextGem_DI_CRISPR_A549_5K"

seqkit_head "$tar_dir/${orig_sample_id}_gex/${orig_sample_id}_gex_S5_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S5_L001_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_gex/${orig_sample_id}_gex_S5_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S5_L001_R2_001.fastq.gz"

seqkit_head "$tar_dir/${orig_sample_id}_crispr/${orig_sample_id}_crispr_S4_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_crispr_subset_S4_L001_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_crispr/${orig_sample_id}_crispr_S4_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_crispr_subset_S4_L001_R2_001.fastq.gz"


# download crispr feature reference
crispr_ref="$raw_dir/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference.csv"
if [[ ! -f "$crisp_ref" ]]; then
  wget "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference.csv" -O "$crispr_ref"
fi

crispr_ref_adjusted="$raw_dir/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference_corrected.csv"
reference_gtf="resources_test/reference_gencodev41_chr1/reference.gtf.gz"
cat "$crispr_ref" | while read line || [[ -n $line ]];
do 
  echo "Line: $line"
  old_id=$( printf "%s\n" "$line" | awk -F',' '{print $7}' )
  echo "Old ID: $old_id"
  if [ "$old_id" = "Non-Targeting" ] || [ "$old_id" = "target_gene_id" ] ; then
    echo "Just writing line"
    printf "%s\n" "$line" >> "$crispr_ref_adjusted"
  else
    gtf_lookup=$(zgrep "$old_id" "$reference_gtf" || test $? = 1;)
    if [ ! -z "$gtf_lookup" ]; then
      echo "Found hit"
      new_id=$(echo "$gtf_lookup" | awk '{if ($3 == "gene") print $10;}' | sed -e "s/^\"//" -e "s/\";$//")
      echo "New ID: $new_id"
      new_line=${line/"$old_id"/"$new_id"}
      printf "%s\n" "$new_line" >> "$crispr_ref_adjusted"
    else
      echo "Did not find hit"
    fi
  fi
done


# Run mapping pipeline
# TODO: Also include conversion to h5mu
cat > /tmp/params.yaml << HERE
param_list:
- id: "$ID"
  input: "$raw_dir"
  library_id:
    - "${orig_sample_id}_gex_subset"
    - "${orig_sample_id}_crispr_subset"
  library_type:
    - "Gene Expression"
    - "CRISPR Guide Capture"

gex_reference: "$genome_tar"
feature_reference: "$crispr_ref_adjusted"
publish_dir: "$OUT/processed"
HERE

nextflow \
  run . \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels.config

# Create h5mu
cat > /tmp/params.yaml << HERE
id: "$ID"
input: "$OUT/processed/10x_5k_lung_crispr.cellranger_multi.output"
publish_dir: "$OUT/"
output: "$orig_sample_id.h5mu"
HERE

nextflow \
  run . \
  -main-script target/nextflow/convert/from_cellranger_multi_to_h5mu/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels.config
