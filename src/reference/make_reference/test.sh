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

echo "> Test succeeded!"
