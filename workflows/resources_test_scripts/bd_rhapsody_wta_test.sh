#!/bin/bash

# settings
ID=bd_rhapsody_wta_test
OUT=resources_test/$ID
DIR="$OUT"
S3DIR=$(echo "$DIR" | sed 's#resources_test#s3://openpipelines-data#')
raw_dir="$OUT/raw"

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# create tempdir
TMPDIR=$(mktemp -d "$VIASH_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# download raw files
mkdir -p "$raw_dir"
wget http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R1_.fastq.gz -O "$raw_dir/sample_R1_.fastq.gz"
wget http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R2_.fastq.gz -O "$raw_dir/sample_R2_.fastq.gz"
wget http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA/GRCh38-PhiX-gencodev29/GRCh38-PhiX-gencodev29-20181205.tar.gz -O "$raw_dir/GRCh38-PhiX-gencodev29-20181205.tar.gz"
wget http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA/GRCh38-PhiX-gencodev29/gencodev29-20181205.gtf -O "$raw_dir/gencodev29-20181205.gtf"

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
  --reference_genome "$raw_dir/GRCh38-PhiX-gencodev29-20181205.tar.gz" \
  --transcriptome_annotation "$raw_dir/gencodev29-20181205.gtf" \
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
