#!/bin/bash



## VIASH START
meta_executable="bin/viash run src/reference/make_reference/config.vsh.yaml --"
## VIASH END

echo "> Running $meta_functionality_name."
fasta="myreference.fa.gz"
gtf="myreference.gtf.gz"

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip

$meta_executable \
  --genome_fasta "GRCh38.primary_assembly.genome.fa.gz" \
  --transcriptome_gtf "gencode.v41.annotation.gtf.gz" \
  --ercc "ERCC92.zip" \
  --subset_regex "(ERCC-00002|chr1)" \
  --output_fasta $fasta \
  --output_gtf $gtf

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f $fasta ]] && echo "Output fasta file could not be found!" && exit 1
[[ ! -f $gtf ]] && echo "Output gtf file could not be found!" && exit 1

echo ">> Checking contents of fasta"
if ! zgrep -q '>chr1' $fasta; then
  echo "Could not find chr1 in output reference!"
  exit 1
fi
if ! zgrep -q '>ERCC-00002' $fasta; then
  echo "Could not find ERCC-00002 in output reference!"
  exit 1
fi
if zgrep -q '>chr11' $fasta; then
  echo "Found chr11 in output reference!"
  exit 1
fi
if zgrep -q '>chr2' $fasta; then
  echo "Found chr2 in output reference!"
  exit 1
fi

echo "> Test succeeded!"