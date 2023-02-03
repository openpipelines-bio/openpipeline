#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# settings
ID=cellranger_tiny_fastq
OUT="resources_test/$ID/"
DIR="$OUT"

# download cellranger tar gz
cellranger_tar_gz="${OUT}/temp_cellranger-6.1.2.tar.gz"
if [ ! -f "$cellranger_tar_gz" ]; then
  echo "Download Cell Ranger 6.1.2 manually first!"
  exit 1
fi

# untar fastqs
cellranger_tiny_fastq="${OUT}/cellranger_tiny_fastq"
if [ ! -f "${cellranger_tiny_fastq}/tinygex_S1_L001_R1_001.fastq.gz" ]; then
  mkdir -p "$cellranger_tiny_fastq"
  
  tar -xzf "$cellranger_tar_gz" \
    -C "$cellranger_tiny_fastq" \
    "cellranger-6.1.2/external/cellranger_tiny_fastq" \
    --strip-components=3
fi

# untar ref
cellranger_tiny_ref="${OUT}/cellranger_tiny_ref"
if [ ! -f "${cellranger_tiny_ref}/reference.json" ]; then
  mkdir -p "$cellranger_tiny_ref"
  
  tar -xzf "$cellranger_tar_gz" \
    -C "$cellranger_tiny_ref" \
    "cellranger-6.1.2/external/cellranger_tiny_ref" \
    --strip-components=3
fi

# Create ref with more recent STAR version
recent_ref_dir="${OUT}/cellranger_tiny_ref_v2_7_10_a"
if [ ! -f "${recent_ref_dir}/Genome" ]; then
  mkdir -p "${recent_ref_dir}"

  target/docker/mapping/star_build_reference/star_build_reference \
    --genome_fasta "$cellranger_tiny_ref/fasta/genome.fa" \
    --output "$recent_ref_dir" \
    --genomeSAindexNbases 7 \
    --transcriptome_gtf "$cellranger_tiny_ref/genes/genes.gtf.gz"
fi

# run cellranger count
bam_dir="${OUT}/bam"
if [ ! -f "$bam_dir/possorted_genome_bam.bam" ]; then
  mkdir -p "$bam_dir"

  viash run src/mapping/cellranger_count/config.vsh.yaml -- \
    --input "$cellranger_tiny_fastq" \
    --reference "$cellranger_tiny_ref" \
    --output "$bam_dir"
fi

# convert to h5mu
raw_h5mu="${OUT}/raw_dataset.h5mu"
if [ ! -f "$step1_h5mu" ]; then
  viash run src/convert/from_10xh5_to_h5mu/config.vsh.yaml -- \
    --input "${bam_dir}/raw_feature_bc_matrix.h5" \
    --output "$raw_h5mu"
fi

# run velocyto
velo_gtf="$cellranger_tiny_ref/genes/genes.gtf.gz"
velo_bam="$bam_dir/possorted_genome_bam.bam"
velo_loom="${OUT}/velocyto.loom"
if [ ! -f "$velo_loom" ]; then
  viash run src/velocity/velocyto/config.vsh.yaml -- \
    --input "$velo_bam" \
    --output "$velo_loom" \
    --transcriptome "$velo_gtf"
fi

# combine raw counts with velocyto data
dataset_h5mu="${OUT}/dataset.h5mu"
if [ ! -f "$dataset_h5mu" ]; then
  viash run src/velocity/velocyto_to_h5mu/config.vsh.yaml -- \
    --input_loom "$velo_loom" \
    --input_h5mu "$raw_h5mu" \
    --output "$dataset_h5mu"
fi

# run htseq
htseq_counts="${OUT}/htseq_counts.tsv"
if [ ! -f "$htseq_counts" ]; then
  viash run src/mapping/htseq_count/config.vsh.yaml -- \
  --input "$velo_bam" \
  --reference "$velo_gtf" \
  --output "$htseq_counts"
fi

multi_star="${OUT}/multi_star"
if [ ! -d "$multi_star" ]; then
  viash run src/mapping/multi_star/config.vsh.yaml -- \
    --input_id "tinygex" \
    --input_r1 "$cellranger_tiny_fastq/tinygex_S1_L001_R1_001.fastq.gz" \
    --input_r2 "$cellranger_tiny_fastq/tinygex_S1_L001_R2_001.fastq.gz" \
    --input_id "tinygex" \
    --input_r1 "$cellranger_tiny_fastq/tinygex_S1_L002_R1_001.fastq.gz" \
    --input_r2 "$cellranger_tiny_fastq/tinygex_S1_L002_R2_001.fastq.gz" \
    --reference_index "$recent_ref_dir" \
    --reference_gtf "$cellranger_tiny_ref/genes/genes.gtf.gz" \
    --output "$multi_star" \
    ---cpus 30
fi