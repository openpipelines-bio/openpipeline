#!/bin/bash

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=bdrhap_5kjrt
OUT=resources_test/$ID

# create raw directory
raw_dir="$OUT/raw"
mkdir -p "$raw_dir"

# Check whether seqkit is available
# TODO: we should turn this into viash components
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

genome_tar="$REPO_ROOT/resources_test/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz"
if [[ ! -f "$genome_tar" ]]; then
    echo "$genome_tar does not exist. Please create the reference genome first"
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

mkdir -p "$raw_dir/GRCh38_primary_assembly_genome_chr1"
cd "$raw_dir/GRCh38_primary_assembly_genome_chr1"
tar -xvf "$genome_tar" 
cd "$REPO_ROOT"

mapping_dir="$raw_dir/mapping_chr_1"
mkdir -p "$mapping_dir"
# MUST USE A STAR THAT IS COMPATIBLE WITH BD RHAPSODY
# For the cwl pipeline 1.9.1, 2.5.2b should work.
docker run --rm -it \
           -v "`pwd`/$OUT:`pwd`/$OUT" \
           -v "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz:$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" \
           -w `pwd` bdgenomics/rhapsody:1.10.1 \
STAR \
    --runThreadN "$n_cores" \
    --genomeDir "$raw_dir/GRCh38_primary_assembly_genome_chr1" \
    --readFilesIn "$tar_dir/12WTA_S1_L432_R1_001.fastq.gz" "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" \
    --runRNGseed 100 \
    --outFileNamePrefix "$mapping_dir/" \
    --readFilesCommand "gzip -d -k -c"

samtools view -F 260 "$mapping_dir/Aligned.out.sam" > "$mapping_dir/primary_aligned_reads.sam"
cut -f 1 "$mapping_dir/primary_aligned_reads.sam" | sort | uniq > "$mapping_dir/mapped_reads.txt"
seqkit grep --threads "$n_cores" -f "$mapping_dir/mapped_reads.txt" "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" > "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq"
seqkit grep --threads "$n_cores" -f "$mapping_dir/mapped_reads.txt" "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" > "$mapping_dir/12WTA_S1_L432_R2_001_chr1.fastq"
gzip -9 -k -c "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq" > "$raw_dir/12WTA_S1_L432_R1_001_chr1.fastq.gz"
gzip -9 -k -c "$mapping_dir/12WTA_S1_L432_R2_001_chr1.fastq" > "$raw_dir/12WTA_S1_L432_R2_001_chr1.fastq.gz"

rm -r "$mapping_dir"

n_reads=100000
seqkit head -n$n_reads "$tar_dir/12SMK_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12SMK_S1_L432_R1_001.fastq.gz"
seqkit head -n$n_reads "$tar_dir/12SMK_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12SMK_S1_L432_R2_001.fastq.gz"
seqkit head -n$n_reads "$tar_dir/12ABC_S1_L432_R1_001.fastq.gz" | gzip -9 > "$raw_dir/12ABC_S1_L432_R1_001.fastq.gz"
seqkit head -n$n_reads "$tar_dir/12ABC_S1_L432_R2_001.fastq.gz" | gzip -9 > "$raw_dir/12ABC_S1_L432_R2_001.fastq.gz"



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



