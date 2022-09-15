#!/bin/bash



## VIASH START
meta_executable="bin/viash run src/reference/make_reference/config.vsh.yaml --"
## VIASH END

echo "> Running $meta_functionality_name."
$meta_executable \
  --genome_fasta "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz" \
  --transcriptome_gtf "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz" \
  --output myreference.tar.gz

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f myreference.fa.gz ]] && echo "Output fasta file could not be found!" && exit 1
[[ ! -f myreference.gtf.gz ]] && echo "Output gtf file could not be found!" && exit 1

echo "> Test succeeded!"