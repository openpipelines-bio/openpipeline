#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=10x_4plex_dtc
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

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}

# dataset page:
# https://www.10xgenomics.com/datasets/40k-mixture-of-dissociated-tumor-cells-from-3-donors-stained-with-totalseqc-antibodies

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/4plex_DTC_kidney_lung_breast_TotalSeqC_fastqs"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"

    # download fastqs and untar
    wget "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/7.2.0/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Multiplex/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Multiplex_fastqs.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

function seqkit_head {
  input="$1"
  output="$2"
  if [[ ! -f "$output" ]]; then
    echo "> Processing `basename $input`"
    seqkit head -n 100000 "$input" | gzip > "$output"
  fi
}

orig_sample_id="4plex_DTC_kidney_lung_breast_TotalSeqC"

seqkit_head "$tar_dir/${orig_sample_id}_gex_S1_L002_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S1_L002_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_gex_S1_L002_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S1_L002_R2_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_gex_S1_L003_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S1_L003_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_gex_S1_L003_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S1_L003_R2_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_gex_S1_L004_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S1_L004_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_gex_S1_L004_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_gex_subset_S1_L004_R2_001.fastq.gz"


seqkit_head "$tar_dir/${orig_sample_id}_ab_S2_L002_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_ab_subset_S2_L002_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_ab_S2_L002_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_ab_subset_S2_L002_R2_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_ab_S2_L003_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_ab_subset_S2_L003_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_ab_S2_L003_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_ab_subset_S2_L003_R2_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_ab_S2_L004_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_ab_subset_S2_L004_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_ab_S2_L004_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_ab_subset_S2_L004_R2_001.fastq.gz"

# download feature reference
feature_ref="$raw_dir/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Multiplex_count_feature_reference.csv"
if [[ ! -f "$feature_ref" ]]; then
  wget "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Multiplex/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Multiplex_count_feature_reference.csv" -O "$feature_ref"
fi

# download probe set
probe_set="$raw_dir/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
if [[ ! -f "$probe_set" ]]; then
  wget "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv" -O "$probe_set"
fi

sed -i 's/#reference_genome=GRCh38/#reference_genome=output/g' "$probe_set"

probe_set_corrected="$raw_dir/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A-corrected.csv"
if [[ ! -f "$probe_set_corrected" ]]; then
  reference_gtf="resources_test/reference_gencodev41_chr1/reference.gtf.gz"
  gunzip -c "$reference_gtf" > "$TMPDIR/uncompressed_ref.gtf" 
  cat "$probe_set" | while read line || [[ -n $line ]];
  do
    echo "Line: $line"
    old_id=$( printf "%s\n" "$line" | awk -F',' '{print $1}' )
    echo "Old ID: $old_id"
    if [[ "$old_id" == "gene_id" ]] || [[ "$old_id" == \#* ]] ; then
      echo "Just writing line"
      printf "%s\n" "$line" >> "$probe_set_corrected"
    else
      gtf_lookup=$(grep "$old_id" "$TMPDIR/uncompressed_ref.gtf" || test $? = 1;)
      if [ ! -z "$gtf_lookup" ]; then
        echo "Found hit"
        new_id=$(echo "$gtf_lookup" | awk '{if ($3 == "gene") print $10;}' | sed -e "s/^\"//" -e "s/\";$//")
        echo "New ID: $new_id"
        new_line=${line/"$old_id"/"$new_id"}
        printf "%s\n" "$new_line" >> "$probe_set_corrected"
      else
        echo "Did not find hit"
      fi
    fi
  done
fi

# Run mapping pipeline
cat > /tmp/params.yaml << HERE
param_list:
- id: "$ID"
  input: "$raw_dir"
  library_id:
    - ${orig_sample_id}_gex_subset
    - ${orig_sample_id}_ab_subset
  library_type:
    - "Gene Expression"
    - "Antibody Capture"
  library_lanes:
    - "any"
    - "any"

probe_set: "$probe_set_corrected"
gex_reference: "$genome_tar"
feature_reference: "$feature_ref"
probe_barcode_ids:
  - BC001+AB001
  - BC002+AB002
  - BC003+AB003
  - BC004+AB004
sample_ids:
  - Kidney_Cancer1_BC1_AB1
  - Kidney_Cancer2_BC2_AB2
  - Breast_Cancer_BC3_AB3
  - Lung_Cancer_BC4_AB4
gex_generate_bam: false
HERE

nextflow \
  run openpipelines-bio/openpipeline \
  -r 3.0.0 \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  --publish_dir "$OUT/processed" \
  -c src/workflows/utils/labels_ci.config

mkdir -p "${OUT}_v10/processed"

nextflow \
  run . \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  --publish_dir "${OUT}_v10/processed" \
  -c src/workflows/utils/labels_ci.config