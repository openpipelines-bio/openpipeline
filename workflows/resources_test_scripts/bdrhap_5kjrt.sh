#!/bin/bash

# TODO: we should turn this into viash components

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=bdrhap_5kjrt
OUT=resources_test/$ID
n_threads=10

# create raw directory
raw_dir="$OUT/raw"
mkdir -p "$raw_dir"

# Check whether seqkit is available
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

# process WTA fastq files
# map to chr1, subsample chr1 reads 


echo "> Processing 12WTA_S1_L432_R[12]_001.fastq.gz"
gzip -d -k -c "$tar_dir/12WTA_S1_L432_R1_001.fastq.gz" > "$raw_dir/12WTA_S1_L432_R1_001.fastq"
gzip -d -k -c "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" > "$raw_dir/12WTA_S1_L432_R2_001.fastq"

mkdir -p "$raw_dir/GRCh38_primary_assembly_genome_chr1"
tar -xvf "$genome_tar" -C "$raw_dir/GRCh38_primary_assembly_genome_chr1"

mapping_dir="$raw_dir/mapping_chr_1"
mkdir -p "$mapping_dir"
# MUST USE A STAR THAT IS COMPATIBLE WITH BD RHAPSODY
# For the cwl pipeline 1.9.1, 2.5.2b should work.
docker run --rm -it \
           -v "`pwd`/$OUT:`pwd`/$OUT" \
           -v "$tar_dir:$tar_dir" \
           -w `pwd` bdgenomics/rhapsody:1.10.1 \
STAR \
    --runThreadN "$n_threads" \
    --genomeDir "$raw_dir/GRCh38_primary_assembly_genome_chr1" \
    --readFilesIn "$tar_dir/12WTA_S1_L432_R1_001.fastq.gz" "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" \
    --runRNGseed 100 \
    --outFileNamePrefix "$mapping_dir/" \
    --readFilesCommand "gzip -d -k -c"

samtools view -F 260 "$mapping_dir/Aligned.out.sam" > "$mapping_dir/primary_aligned_reads.sam"
cut -f 1 "$mapping_dir/primary_aligned_reads.sam" | sort | uniq > "$mapping_dir/mapped_reads.txt"
seqkit grep --threads "$n_threads" -f "$mapping_dir/mapped_reads.txt" "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" > "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq"
seqkit grep --threads "$n_threads" -f "$mapping_dir/mapped_reads.txt" "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" > "$mapping_dir/12WTA_S1_L432_R2_001_chr1.fastq"
gzip -9 -k -c "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq" > "$raw_dir/12WTA_S1_L432_R1_001_chr1.fastq.gz"
gzip -9 -k -c "$mapping_dir/12WTA_S1_L432_R2_001_chr1.fastq" > "$raw_dir/12WTA_S1_L432_R2_001_chr1.fastq.gz"

rm -r "$mapping_dir"

# subsample other files
echo "> Processing 12SMK_S1_L432_R1_001.fastq.gz"
seqkit head -n 5000000 "$tar_dir/12SMK_S1_L432_R1_001.fastq.gz" | gzip > "$raw_dir/12SMK_S1_L432_R1_001.fastq.gz"
echo "> Processing 12SMK_S1_L432_R2_001.fastq.gz"
seqkit head -n 5000000 "$tar_dir/12SMK_S1_L432_R2_001.fastq.gz" | gzip > "$raw_dir/12SMK_S1_L432_R2_001.fastq.gz"
echo "> Processing 12ABC_S1_L432_R1_001.fastq.gz"
seqkit head -n 5000000 "$tar_dir/12ABC_S1_L432_R1_001.fastq.gz" | gzip > "$raw_dir/12ABC_S1_L432_R1_001.fastq.gz"
echo "> Processing 12ABC_S1_L432_R2_001.fastq.gz"
seqkit head -n 5000000 "$tar_dir/12ABC_S1_L432_R2_001.fastq.gz" | gzip > "$raw_dir/12ABC_S1_L432_R2_001.fastq.gz"

# copy immune panel fasta
cp "$tar_dir/BDAbSeq_ImmuneDiscoveryPanel.fasta" "$raw_dir"



# process samples with bd rhap component
# TODO: change to bd rhap ingestion pipeline
bdrhap_5kjrt="resources_test/bdrhap_5kjrt/raw"
bdrhap_ref_gencodev40_chr1="resources_test/bdrhap_ref_gencodev40_chr1"
cat > /tmp/params.yaml << HERE
param_list:
- id: "SMK"
  run_name: "SMK"
  input: "$bdrhap_5kjrt/12SMK_S1_L432_R[12]_001.fastq.gz"
  sample_tags_version: "hs"
  tag_names: ["1-Jurkat", "2-Ramos", "3-THP1"]
- id: "ABC"
  run_name: "ABC"
  input: "$bdrhap_5kjrt/12ABC_S1_L432_R[12]_001.fastq.gz"
  abseq_reference: "$bdrhap_5kjrt/BDAbSeq_ImmuneDiscoveryPanel.fasta"
- id: "WTA"
  run_name: "WTA"
  input: "$bdrhap_5kjrt/12WTA_S1_L432_R[12]_001.fastq.gz"
reference_genome: "$bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz"
transcriptome_annotation: "$bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf"
publish_dir: "output_test"
putative_cell_call: "mRNA"
exact_cell_count: 500
HERE

bin/nextflow \
  run . \
  -main-script target/nextflow/mapping/bd_rhapsody_wta/main.nf \
  -resume \
  -profile docker,mount_temp \
  -with-trace work/trace.txt \
  -params-file /tmp/params.yaml \
  -c workflows/utils/labels.config \
  -c workflows/utils/errorstrat_ignore.config
