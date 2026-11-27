#!/bin/bash

set -eo pipefail

## VIASH START
meta_executable="bin/viash run src/reference/make_reference/config.vsh.yaml --"
## VIASH END

# Fetch test data
echo ">> Fetching test data"

wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz
wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip

# Check that the GTF line of an ERCC-00002 feature contains an expected attribute
# Usage: assert_ercc_attribute <gtf.gz> <feature> <attribute>
function assert_ercc_attribute {
  local gtf="$1" feature="$2" attribute="$3"
  local line
  line=$(zcat "$gtf" | awk -F'\t' -v feature="$feature" '$1 == "ERCC-00002" && $3 == feature')
  if [[ -z "$line" ]]; then
    echo "Could not find a '$feature' feature for ERCC-00002 in $gtf!"
    exit 1
  fi
  if [[ "$line" != *"$attribute"* ]]; then
    echo "The '$feature' feature for ERCC-00002 does not contain '$attribute':"
    echo "$line"
    exit 1
  fi
}

# Test 1
echo ">> Test1"
mkdir test1
pushd test1
fasta="myreference.fa.gz"
gtf="myreference.gtf.gz"

"$meta_executable" \
  --genome_fasta "../Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz" \
  --transcriptome_gtf "../Homo_sapiens.GRCh38.109.chr.gtf.gz" \
  --ercc "../ERCC92.zip" \
  --subset_regex "(ERCC-00002|1)" \
  --output_fasta $fasta \
  --output_gtf $gtf

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f $fasta ]] && echo "Output fasta file could not be found!" && exit 1
[[ ! -f $gtf ]] && echo "Output gtf file could not be found!" && exit 1

echo ">> Checking contents of fasta"
if ! zgrep -q '>1' $fasta; then
  echo "Could not find chromosome '1' in output reference!"
  exit 1
fi
if ! zgrep -q '>ERCC-00002' $fasta; then
  echo "Could not find ERCC-00002 in output reference!"
  exit 1
fi

echo ">> Checking ERCC features and attributes in gtf"
for feature in gene transcript exon; do
  assert_ercc_attribute $gtf $feature 'gene_name "ERCC-00002";'
  assert_ercc_attribute $gtf $feature 'gene_type "spike_in_sequence";'
done
for feature in transcript exon; do
  assert_ercc_attribute $gtf $feature 'transcript_type "spike_in_sequence";'
done
popd

# Test 2
echo ">> Test 2"
mkdir test2
pushd test2
fasta="myreference.fa.gz"
gtf="myreference.gtf.gz"

"$meta_executable" \
  --genome_fasta "../Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz" \
  --transcriptome_gtf "../Homo_sapiens.GRCh38.109.chr.gtf.gz" \
  --output_fasta $fasta \
  --output_gtf $gtf

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f $fasta ]] && echo "Output fasta file could not be found!" && exit 1
[[ ! -f $gtf ]] && echo "Output gtf file could not be found!" && exit 1

echo ">> Checking contents of fasta"
if ! zgrep -q '>1' $fasta; then
  echo "Could not find chromosome '1' in output reference!"
  exit 1
fi
if zgrep -q '>ERCC-00002' $fasta; then
  echo "Should not find ERCC-00002 in output reference!"
  exit 1
fi
popd

# Small inputs for the ERCC annotation tests
printf '>1\nACGTACGTACGTACGTACGT\n' > small_genome.fa
printf '1\ttest\tgene\t1\t20\t.\t+\t.\tgene_id "G1"; gene_name "G1"; gene_type "protein_coding";\n' > small_annotation.gtf

# ERCC archive in which ERCC-00002 already has gene_type and transcript_type attributes
mkdir ercc_typed
unzip -q ERCC92.zip -d ercc_typed
sed -i '/^ERCC-00002\t/ s/$/ gene_type "custom_gene"; transcript_type "custom_transcript";/' ercc_typed/ERCC92.gtf
(cd ercc_typed && python3 -m zipfile -c ../ERCC92_typed.zip ERCC92.fa ERCC92.gtf)

# Test 3
echo ">> Test 3: custom ERCC gene and transcript types"
mkdir test3
pushd test3
gtf="myreference.gtf.gz"

"$meta_executable" \
  --genome_fasta "../small_genome.fa" \
  --transcriptome_gtf "../small_annotation.gtf" \
  --ercc "../ERCC92.zip" \
  --ercc_gene_type "my_gene_type" \
  --ercc_transcript_type "my_transcript_type" \
  --output_fasta myreference.fa.gz \
  --output_gtf $gtf

echo ">> Checking ERCC types in gtf"
for feature in gene transcript exon; do
  assert_ercc_attribute $gtf $feature 'gene_type "my_gene_type";'
done
for feature in transcript exon; do
  assert_ercc_attribute $gtf $feature 'transcript_type "my_transcript_type";'
done
popd

# Test 4
echo ">> Test 4: existing ERCC types are kept by default"
mkdir test4
pushd test4
gtf="myreference.gtf.gz"

"$meta_executable" \
  --genome_fasta "../small_genome.fa" \
  --transcriptome_gtf "../small_annotation.gtf" \
  --ercc "../ERCC92_typed.zip" \
  --output_fasta myreference.fa.gz \
  --output_gtf $gtf

echo ">> Checking ERCC types in gtf"
for feature in gene transcript exon; do
  assert_ercc_attribute $gtf $feature 'gene_type "custom_gene";'
done
for feature in transcript exon; do
  assert_ercc_attribute $gtf $feature 'transcript_type "custom_transcript";'
done
popd

# Test 5
echo ">> Test 5: existing ERCC types are replaced with --ercc_overwrite_types"
mkdir test5
pushd test5
gtf="myreference.gtf.gz"

"$meta_executable" \
  --genome_fasta "../small_genome.fa" \
  --transcriptome_gtf "../small_annotation.gtf" \
  --ercc "../ERCC92_typed.zip" \
  --ercc_overwrite_types \
  --output_fasta myreference.fa.gz \
  --output_gtf $gtf

echo ">> Checking ERCC types in gtf"
for feature in gene transcript exon; do
  assert_ercc_attribute $gtf $feature 'gene_type "spike_in_sequence";'
done
for feature in transcript exon; do
  assert_ercc_attribute $gtf $feature 'transcript_type "spike_in_sequence";'
done
popd

echo "> Test succeeded!"
