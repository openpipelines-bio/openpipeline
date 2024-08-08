#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=bdrhap_5kjrt
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

# check whether reference is available
reference_dir="resources_test/reference_gencodev41_chr1"
genome_tar="$reference_dir/reference_bd_rhapsody.tar.gz"
if [[ ! -f "$genome_tar" ]]; then
    echo "$genome_tar does not exist. Please create the reference genome first"
    exit 1
fi

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/12WTA-ABC-SMK-EB-5kJRT"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"
    wget "http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-Demo-Data-Inputs/12WTA-ABC-SMK-EB-5kJRT.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

genome_dir="$raw_dir/temp_reference_gencodev41_chr1"
if [[ ! -d "$genome_dir" ]]; then
  echo "> Untarring genome"
  mkdir -p "$genome_dir"
  tar -xvf "$genome_tar" -C "$genome_dir"
fi

# process WTA fastq files
# map to chr1, subsample chr1 reads 
mapping_dir="$raw_dir/temp_mapping_chr_1"
if [[ ! -f "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq" ]]; then
  echo "> Processing 12WTA_S1_L432_R[12]_001.fastq.gz"
  mkdir -p "$mapping_dir"
  # MUST USE A STAR THAT IS COMPATIBLE WITH BD RHAPSODY
  # For the cwl pipeline 1.9.1, 2.5.2b should work.
  echo "star"
  docker run --rm -i \
    -v "`pwd`/$OUT:`pwd`/$OUT" \
    -v "$tar_dir:$tar_dir" \
    -w `pwd` bdgenomics/rhapsody:1.10.1 \
    STAR \
      --runThreadN "$n_threads" \
      --genomeDir "$genome_dir" \
      --readFilesIn "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" \
      --runRNGseed 100 \
      --outFileNamePrefix "$mapping_dir/" \
      --readFilesCommand "gzip -d -k -c" \
      --clip3pAdapterSeq "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" \
      --outFilterMatchNmin "25" \
      --quantTranscriptomeBan "Singleend" # Prohibit mapping of one side of the read
  # chown to current user before removing mapping dir
  docker run --rm -i -v "`pwd`/$OUT:`pwd`/$OUT" -w `pwd` bdgenomics/rhapsody:1.10.1 \
    chown "$(id -u):$(id -g)" --silent --recursive "$mapping_dir/"

  echo "samtools"
  samtools view -F 260 "$mapping_dir/Aligned.out.sam" > "$mapping_dir/primary_aligned_reads.sam"
  echo "cut"
  cut -f 1 "$mapping_dir/primary_aligned_reads.sam" | sort | uniq > "$mapping_dir/mapped_reads.txt"
  head -500000 "$mapping_dir/mapped_reads.txt" > "$mapping_dir/mapped_reads_subset.txt"
  echo "seqkit"
  seqkit grep --threads "$n_threads" -f "$mapping_dir/mapped_reads_subset.txt" "$tar_dir/12WTA_S1_L432_R1_001.fastq.gz" > "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq"
  seqkit grep --threads "$n_threads" -f "$mapping_dir/mapped_reads_subset.txt" "$tar_dir/12WTA_S1_L432_R2_001.fastq.gz" > "$mapping_dir/12WTA_S1_L432_R2_001_chr1.fastq"

  # rm -r "$mapping_dir"
  # rm -r "$genome_dir"
fi

# subsample other files
smk_r1_file="$raw_dir/12SMK_S1_L432_R1_001.fastq.gz"
if [[ ! -f "$smk_r1_file" ]]; then
  echo "> Processing `basename $smk_r1_file`"
  cp "$tar_dir/12SMK_S1_L432_R1_001.fastq.gz" "$smk_r1_file"
fi
smk_r2_file="$raw_dir/12SMK_S1_L432_R2_001.fastq.gz"
if [[ ! -f "$smk_r2_file" ]]; then
  echo "> Processing `basename $smk_r2_file`"
  cp "$tar_dir/12SMK_S1_L432_R2_001.fastq.gz" "$smk_r2_file"
fi
abc_r1_file="$raw_dir/12ABC_S1_L432_R1_001_subset.fastq.gz"
if [[ ! -f "$abc_r1_file" ]]; then
  echo "> Processing `basename $abc_r1_file`"
  seqkit head -n 500000 "$tar_dir/12ABC_S1_L432_R1_001.fastq.gz" | gzip > "$abc_r1_file"
fi
abc_r2_file="$raw_dir/12ABC_S1_L432_R2_001_subset.fastq.gz"
if [[ ! -f "$abc_r2_file" ]]; then
  echo "> Processing `basename $abc_r2_file`"
  seqkit head -n 500000 "$tar_dir/12ABC_S1_L432_R2_001.fastq.gz" | gzip > "$abc_r2_file"
fi
wta_r1_file="$raw_dir/12WTA_S1_L432_R1_001_subset.fastq.gz"
if [[ ! -f "$wta_r1_file" ]]; then
  echo "> Processing `basename $wta_r1_file`"
  gzip -9 -k -c "$mapping_dir/12WTA_S1_L432_R1_001_chr1.fastq" > "$wta_r1_file"
fi
wta_r2_file="$raw_dir/12WTA_S1_L432_R2_001_subset.fastq.gz"
if [[ ! -f "$wta_r2_file" ]]; then
  echo "> Processing `basename $wta_r2_file`"
  gzip -9 -k -c "$mapping_dir/12WTA_S1_L432_R2_001_chr1.fastq" > "$wta_r2_file"
fi
# copy immune panel fasta
fasta_file="$raw_dir/BDAbSeq_ImmuneDiscoveryPanel.fasta"
if [[ ! -f "$fasta_file" ]]; then
  cp "$tar_dir/BDAbSeq_ImmuneDiscoveryPanel.fasta" "$fasta_file"
fi


# process samples with bd rhap component
# TODO: change to bd rhap ingestion pipeline
cat > /tmp/params.yaml << HERE
param_list:
- id: "SMK"
  input: "$smk_r1_file;$smk_r2_file"
  sample_tags_version: "hs"
  tag_names: ["1-Jurkat", "2-Ramos", "3-THP1"]
- id: "ABC"
  input: "$abc_r1_file;$abc_r2_file"
  abseq_reference: "$fasta_file"
- id: "WTA"
  input: "$wta_r1_file;$wta_r2_file"
mode: wta
reference: "$reference_dir/reference_bd_rhapsody.tar.gz"
transcriptome_annotation: "$reference_dir/reference.gtf.gz"
publish_dir: "$OUT/processed"
putative_cell_call: "mRNA"
exact_cell_count: 4000
HERE

nextflow \
  run . \
  -main-script src/workflows/ingestion/bd_rhapsody/main.nf \
  -resume \
  -profile docker,mount_temp \
  -with-trace work/trace.txt \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels.config \
  -c src/workflows/utils/errorstrat_ignore.config


wta_reads= "$raw_dir/12WTA_S1_L432_R1_001_subset.fastq.gz;$raw_dir/12WTA_S1_L432_R2_001_subset.fastq.gz"
abc_reads= "$raw_dir/12ABC_S1_L432_R1_001_subset.fastq.gz;$raw_dir/12ABC_S1_L432_R2_001_subset.fastq.gz"
smk_reads= "$raw_dir/12SMK_S1_L432_R1_001.fastq.gz;$raw_dir/12SMK_S1_L432_R2_001.fastq.gz"

nextflow \
  run . \
  -main-script target/nextflow/workflows/ingestion/bd_rhapsody2/main.nf  \
  -resume \
  -profile docker,mount_temp \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/errorstrat_ignore.config \
  --reads "$wta_reads;$abc_reads;$smk_reads" \
  --reference_archive "$reference_dir/reference_bd_rhapsody_v2.tar.gz" \
  --abseq_reference "$raw_dir/BDAbSeq_ImmuneDiscoveryPanel.fasta" \
  --output_dir "processed2" \
  --cell_calling_data "mRNA" \
  --exact_cell_count 4900 \
  --publish_dir "$OUT"
