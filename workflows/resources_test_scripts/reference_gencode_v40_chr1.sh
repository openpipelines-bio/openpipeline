#!/bin/bash

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=reference_gencode_v40_chr1
OUT=resources_test/$ID

# Check whether STAR and seqkit are available
# TODO: we should turn this into viash components
if ! command -v STAR &> /dev/null; then
    echo "This script requires STAR. Please make sure the binary is added to your PATH."
    exit 1
fi

if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

# create output dir
mkdir -p "$OUT"

## gtf
# download raw files and subset reference to chromosome 1
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz" -O "$OUT/gencode_v40_annotation.gtf.gz" 
gunzip -c "$OUT/gencode_v40_annotation.gtf.gz" > "$OUT/gencode_v40_annotation.gtf"
grep -E '^(##|chr1[^0-9])' "$OUT/gencode_v40_annotation.gtf" > "$OUT/gencode_v40_annotation_chr1.gtf"

# remove temp files
rm "$OUT/gencode_v40_annotation.gtf.gz"
rm "$OUT/gencode_v40_annotation.gtf"

## fasta
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz" -O "$OUT/GRCh38_primary_assembly_genome.fa.gz"
gunzip -c "$OUT/GRCh38_primary_assembly_genome.fa.gz" > "$OUT/GRCh38_primary_assembly_genome.fa"
seqkit grep -r -p 'chr1 1' -n "$OUT/GRCh38_primary_assembly_genome.fa" > "$OUT/GRCh38_primary_assembly_genome_chr1.fa"
gzip -9 -k "$OUT/GRCh38_primary_assembly_genome_chr1.fa"

# run STAR to generate reference compatible with BD rhapsody
# MUST USE A STAR THAT IS COMPATIBLE WITH BD RHAPSODY
# For the cwl pipeline 1.9.1, 2.5.2b should work.
mkdir "$OUT/GRCh38_primary_assembly_genome_chr1"
STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir "$OUT/GRCh38_primary_assembly_genome_chr1" \
  --genomeFastaFiles "$OUT/GRCh38_primary_assembly_genome_chr1.fa" \
  --sjdbGTFfile "$OUT/gencode_v40_annotation.gtf" \
  --sjdbOverhang 100 \
  --genomeSAindexNbases 11
tar -czf "$OUT/GRCh38_primary_assembly_genome_chr1.tar.gz" "$OUT/GRCh38_primary_assembly_genome_chr1"

rm -r "$OUT/GRCh38_primary_assembly_genome_chr1"
rm "$OUT/GRCh38_primary_assembly_genome_chr1.fa.gz"
rm "$OUT/GRCh38_primary_assembly_genome_chr1.fa"
rm "$OUT/GRCh38_primary_assembly_genome.fa.gz"
rm "$OUT/GRCh38_primary_assembly_genome.fa"
