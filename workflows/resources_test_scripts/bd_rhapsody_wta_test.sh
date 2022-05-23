#!/bin/bash

# settings
ID=bd_rhapsody_wta_test
OUT=resources_test/$ID
DIR="$OUT"
S3DIR=$(echo "$DIR" | sed 's#resources_test#s3://openpipelines-data#')
raw_dir="$OUT/raw"

# Check if STAR and seqkit are available
if ! command -v STAR &> /dev/null
then
    echo "This script requires STAR. Please make sure the binary is added to your PATH."
    exit
fi

if ! command -v seqkit &> /dev/null
then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit
fi

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# create tempdir
TMPDIR=$(mktemp -d "$VIASH_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# download raw files and subset reference to chromsosome 1
mkdir -p "$raw_dir"
wget http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R1_.fastq.gz -O "$raw_dir/sample_R1_.fastq.gz"
wget http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R2_.fastq.gz -O "$raw_dir/sample_R2_.fastq.gz"
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz -O "$raw_dir/gencode_v40_annotation.gtf.gz" 
gunzip -c "$raw_dir/gencode_v40_annotation.gtf.gz" > "$raw_dir/gencode_v40_annotation.gtf" && rm "$raw_dir/gencode_v40_annotation.gtf.gz"
grep  -E '^(##|chr1\t)' "$raw_dir/gencode_v40_annotation.gtf" > "$raw_dir/gencode_v40_annotation_chr1.gtf"

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz -O "$raw_dir/GRCh38_primary_assembly_genome.fa.gz"
gunzip -c "$raw_dir/GRCh38_primary_assembly_genome.fa.gz" > "$raw_dir/GRCh38_primary_assembly_genome.fa"
seqkit grep -r -p 'chr1 1' -n "$raw_dir/GRCh38_primary_assembly_genome.fa" > "$raw_dir/GRCh38_primary_assembly_genome_chr1.fa"
gzip -k "$raw_dir/GRCh38_primary_assembly_genome_chr1.fa"

# run STAR to generate reference compatible with BD rhapsody
# MUST USE A STAR THAT IS COMPATIBLE WITH BD RHAPSODY
# For the cwl pipeline 1.9.1, 2.5.2b should work.
# If you get an error regarding resources, manually change all values in the '"requirements":' field to 
# values within the bound of your machine's CPU en memory specifications.  
mkdir "$raw_dir/GRCh38_primary_assembly_genome_chr1"
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir "$raw_dir/GRCh38_primary_assembly_genome_chr1" --genomeFastaFiles "$raw_dir/GRCh38_primary_assembly_genome_chr1.fa" --sjdbGTFfile "$raw_dir/gencode_v40_annotation.gtf" --sjdbOverhang 100 --genomeSAindexNbases 11
tar -czf "$raw_dir/GRCh38_primary_assembly_genome_chr1.tar.gz" "$raw_dir/GRCh38_primary_assembly_genome_chr1" && rm -r "$raw_dir/GRCh38_primary_assembly_genome_chr1"

rm "$raw_dir/GRCh38_primary_assembly_genome.fa.gz"
rm "$raw_dir/GRCh38_primary_assembly_genome.fa"

# create csvs
cat > "$raw_dir/input.csv" << HERE
id,input
sample_RSEC,sample_R*_.fastq.gz
HERE

cat > "$raw_dir/input_remote.csv" << HERE
id,input
sample_RSEC,http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R1_.fastq.gz;http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R2_.fastq.gz
HERE

# process raw files
processed_dir="$OUT/processed"
mkdir -p "$processed_dir"

# run bdrhap
bdrhap_out="$processed_dir/bdrhap_out/"
mkdir -p "$bdrhap_out"

target/docker/mapping/bd_rhapsody_wta/bd_rhapsody_wta \
  --input "$raw_dir/sample_R1_.fastq.gz" \
  --input "$raw_dir/sample_R2_.fastq.gz" \
  --reference_genome "$raw_dir/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  --transcriptome_annotation "$raw_dir/gencode_v40_annotation_chr1.gtf" \
  --output "$bdrhap_out"

# convert to h5ad
h5ad_out="$processed_dir/output.h5ad"
target/docker/convert/from_bdrhap_to_h5ad/from_bdrhap_to_h5ad \
  --input "$bdrhap_out" \
  --id "sample_RSEC" \
  --output "$h5ad_out"
  src/convert/from_bdrhap_to_h5ad/config.vsh.yaml

# nextflow \
#   run . \
#   -main-script workflows/1_ingestion/bd_rhapsody_wta/main.nf \
#   --csv "$raw_dir/input.csv" \
#   --reference_genome "$raw_dir/GRCh38-PhiX-gencodev29-20181205.tar.gz" \
#   --transcriptome_annotation "$raw_dir/gencodev29-20181205.gtf" \
#   --output "$processed_dir" \
#   -resume


# aws s3 sync --profile xxx "$DIR" "$S3DIR"
