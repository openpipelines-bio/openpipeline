#!/bin/bash

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=bdrhap_ref_gencodev41_chr1
OUT=resources_test/$ID
n_threads=30

# get reference
viash run src/reference/make_reference/config.vsh.yaml -- \
  --genome_fasta https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz \
  --transcriptome_gtf https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz \
  --subset_regex chr1 \
  --output_fasta "$OUT/reference_gencode_v41_chr1.fa.gz" \
  --output_gtf "$OUT/reference_gencode_v41_chr1.gtf.gz"

# build as bdrhap star index
viash run src/reference/build_bdrhap_reference/config.vsh.yaml -- \
  --genome_fasta "$OUT/reference_gencode_v41_chr1.fa.gz" \
  --transcriptome_gtf "$OUT/reference_gencode_v41_chr1.gtf.gz" \
  --output "$OUT/gencode_v41_annotation_star.tar.gz"
