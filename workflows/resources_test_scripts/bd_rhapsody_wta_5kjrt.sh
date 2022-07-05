#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=bd_rhapsody_wta_5kjrt
OUT=resources_test/$ID
raw_dir="$OUT/raw"

# create output directory
mkdir -p "$raw_dir"

# Check whether seqkit is available
# TODO: we should turn this into viash components
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

if [[ ! -f "resources_test/reference_gencode_v40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz" ]]; then
    echo "resources_test/reference_gencode_v40_chr1/GRCh38_primary_assembly_genome_chr1 does not exist. Please create the reference genome first"
    exit 1
fi

tar_dir="$HOME/.cache/openpipeline/12WTA-ABC-SMK-EB-5kJRT"

if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"
    wget "http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-Demo-Data-Inputs/12WTA-ABC-SMK-EB-5kJRT.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

# process files 
n_cores=10

gzip -d -k -c "$tar_dir/12WTA_S1_L432_R1_001.fastq.gz" > "$raw_dir/12WTA_S1_L432_R1_001.fastq"
gzip -d -k -c "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" > "$raw_dir/12WTA_S1_L432_R2_001.fastq"

mkdir "$raw_dir/GRCh38_primary_assembly_genome_chr1"
cd "$raw_dir/GRCh38_primary_assembly_genome_chr1"
tar -xvf "$REPO_ROOT/resources_test/reference_gencode_v40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz" 
cd "$REPO_ROOT"

mapping_dir="$raw_dir/mapping_chr_1"
mkdir -p "$mapping_dir"
STAR \
    --runThreadN 8 \
    --genomeDir "resources_test/bd_rhapsody_wta_5kjrt/raw/GRCh38_primary_assembly_genome_chr1" \
    --readFilesIn "$raw_dir/12WTA_S1_L432_R1_001.fastq" "$raw_dir/12WTA_S1_L432_R2_001.fastq" \
    --runRNGseed 100 \
    --outFileNamePrefix "$mapping_dir"

samtools view -F 260 "$mapping_dir/Aligned.out.sam" > "$mapping_dir/primary_aligned_reads.sam"
cut -f 1 "$mapping_dir/primary_aligned_reads.sam" | sort | uniq > "$mapping_dir/mapped_reads.txt"
seqkit grep -f "$mapping_dir/mapped_reads.txt" "$mapping_dir/12WTA_S1_L432_R1_001.fastq" > "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq"
seqkit grep -f "$mapping_dir/mapped_reads.txt" "$raw_dir/12WTA_S1_L432_R2_001.fastq" > "$mapping_dir/12WTA_S1_L432_R2_001_chr1.fastq"
gzip -9 -k -c "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq" > "$raw_dir/12WTA_S1_L432_R1_001_chr1.fastq.gz"
gzip -9 -k -c "$mapping_dir/12WTA_S1_L432_R2_001_chr1.fastq" > "$raw_dir/12WTA_S1_L432_R2_001_chr1.fastq.gz"

rm -r "$mapping_dir"

seqkit head -n100 "$tar_dir/12SMK_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12SMK_S1_L432_R1_001.fastq.gz"
seqkit head -n100 "$tar_dir/12SMK_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12SMK_S1_L432_R2_001.fastq.gz"
seqkit head -n100 "$tar_dir/12ABC_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12ABC_S1_L432_R1_001.fastq.gz"
seqkit head -n100 "$tar_dir/12ABC_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12ABC_S1_L432_R2_001.fastq.gz"


cp "$tar_dir/BDAbSeq_ImmuneDiscoveryPanel.fasta" "$raw_dir"


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



