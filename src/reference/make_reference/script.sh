#!/bin/bash

set -eo pipefail

## VIASH START
par_genome_fasta="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz"
par_transcriptome_gtf="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz"
par_ercc="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip"
# par_subset_regex='(ERCC-00002|chr1)'
par_output_fasta="temp/reference.fa.gz"
par_output_gtf="temp/reference.gtf.gz"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

echo "> Getting path of fasta file"
par_genome_fasta=$(realpath $par_genome_fasta)
echo "> Getting path of annotation file"
par_transcriptome_gtf=$(realpath $par_transcriptome_gtf)

echo "> Processing genome sequence"
genome_fasta="$tmpdir/genome_sequence.fa"
# if genome is gzipped, extract. otherwise not
if file --mime-type "$par_genome_fasta" | grep -q gzip$; then
  zcat "$par_genome_fasta" > "$genome_fasta"
else
  cp "$par_genome_fasta" "$genome_fasta"
fi

echo "> Processing transcriptome annotation"
transcriptome_gtf="$tmpdir/transcriptome_annotation.gtf"
# if transcriptome is gzipped, extract. otherwise not
if file --mime-type "$par_transcriptome_gtf" | grep -q gzip$; then
  zcat "$par_transcriptome_gtf" > "$transcriptome_gtf"
else
  cp "$par_transcriptome_gtf" "$transcriptome_gtf"
fi

if [[ ! -z $par_ercc ]]; then
  echo "> Processing ERCC sequences"
  # wget "$par_ercc" -O "$tmpdir/ercc.zip"
  # unzip "$tmpdir/ercc.zip" -d "$tmpdir"
  unzip "$par_ercc" -d "$tmpdir"
  cat "$tmpdir/ERCC92.fa" >> "$genome_fasta"
  cat "$tmpdir/ERCC92.gtf" >> "$transcriptome_gtf"
fi

# create output & filter reference if so desired
if [[ ! -z $par_subset_regex ]]; then
  echo "> Subsetting reference with regex '$par_subset_regex'"
  awk '{print $1}' "$genome_fasta" | seqkit grep -r -p "^$par_subset_regex\$" > "$tmpdir/genome_sequence_filtered.fa"
  genome_fasta="$tmpdir/genome_sequence_filtered.fa"
  grep -E "^$par_subset_regex[^A-Za-z0-9]" "$transcriptome_gtf" > "$tmpdir/transcriptome_annotation_filtered.gtf"
  transcriptome_gtf="$tmpdir/transcriptome_annotation_filtered.gtf"

  echo
  echo "Matched tags:"
  cat "$genome_fasta" | grep '^>' | sed 's#^>##' | sed 's# .*##' | sort | uniq
  echo 
fi

echo "> Gzipping outputs"
pigz -c "$genome_fasta" > "$par_output_fasta"
pigz -c "$transcriptome_gtf" > "$par_output_gtf"

# to do: re enable
# echo "> Sanity check of outputs"
# readarray -t fasta_tags < <( cat "$genome_fasta" | grep '^>' | sed 's#^>##' | sed 's# .*##' | sort | uniq )
# readarray -t transcriptome_tags < <( cat "$transcriptome_gtf" | cut -d$'\t' -f1 | sort | uniq | grep '^[^#]' )
# [ "${fasta_tags[*]}" == "${transcriptome_tags[*]}" ] || { echo "Warning: fasta tags differ from transcriptome tags"; exit 1; }