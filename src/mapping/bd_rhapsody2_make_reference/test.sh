#!/bin/bash

set -e

############################################
# download test data
TMP_DIR=/tmp/bd_rhapsody_make_reference

ORIG_FA=$TMP_DIR/reference.fa.gz
ORIG_GTF=$TMP_DIR/reference.gtf.gz


# create temporary directory and clean up on exit
mkdir -p $TMP_DIR
function clean_up {
    rm -rf "$TMP_DIR"
}
trap clean_up EXIT

# fetch reference
if [ ! -f $ORIG_FA ]; then
  wget --no-check-certificate https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz \
    -O $ORIG_FA
fi
if [ ! -f $ORIG_GTF ]; then
  wget --no-check-certificate https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz \
    -O $ORIG_GTF
fi

# create small reference
START=30000
END=31500
CHR=chr1

# subset to small region
seqkit grep -r -p "^$CHR\$" "$ORIG_FA" | \
  seqkit subseq -r "$START:$END" > $OUT_DIR/reference_small.fa

zcat "$ORIG_GTF" | \
  awk -v FS='\t' -v OFS='\t' "
    \$1 == \"$CHR\" && \$4 >= $START && \$5 <= $END {
      \$4 = \$4 - $START + 1;
      \$5 = \$5 - $START + 1;
      print;
    }" > $OUT_DIR/reference_small.gtf

############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_doesnt_exist() {
  [ ! -f "$1" ] || { echo "File '$1' exists but shouldn't" && exit 1; }
}
assert_file_empty() {
  #  () will execute in a shubshell, could you use {;}?
  [ ! -s "$1" ] || { echo "File '$1' is not empty but should be" && exit 1; }
}
assert_file_not_empty() {
  # [ -s "$1" ] || (echo "File '$1' is empty but shouldn't be" && exit 1)
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  # grep -q "$2" "$1" || (echo "File '$1' does not contain '$2'" && exit 1)
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains() {
  # grep -q "$2" "$1" && (echo "File '$1' contains '$2' but shouldn't" && exit 1)
  grep -q "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}

echo "#############################################"
echo "> Simple run"

mkdir simple_run
cd simple_run

out_tar="myreference.tar.gz"

echo "> Running $meta_name."
$meta_executable \
  --genome_fasta "$ORIG_FA" \
  --gtf "$ORIG_GTF" \
  --reference_archive "$out_tar" \
  --extra_star_params "--genomeSAindexNbases 6" \
  ---cpus 2

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

assert_file_exists "$out_tar"
assert_file_not_empty "$out_tar"

echo ">> Checking whether output contains the expected files"
tar -xvf "$out_tar" > /dev/null
assert_file_exists "BD_Rhapsody_Reference_Files/star_index/genomeParameters.txt"
assert_file_exists "BD_Rhapsody_Reference_Files/bwa-mem2_index/reference_small.ann"
assert_file_exists "BD_Rhapsody_Reference_Files/reference_small-processed.gtf"
assert_file_exists "BD_Rhapsody_Reference_Files/mitochondrial_contigs.txt"
assert_file_contains "BD_Rhapsody_Reference_Files/reference_small-processed.gtf" "chr1.*HAVANA.*ENSG00000243485"
assert_file_contains "BD_Rhapsody_Reference_Files/mitochondrial_contigs.txt" 'chrMT'

cd ..

echo "#############################################"

echo "> Tests succeeded!"