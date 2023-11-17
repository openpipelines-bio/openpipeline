#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=reference_gencodev41_chr1
OUT=resources_test/$ID


wget "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip" -O "$OUT/ERCC92.zip"

NXF_VER=21.10.6 nextflow \
  run . \
  -main-script src/workflows/ingestion/make_reference/main.nf \
  -profile docker \
  --id "$ID" \
  --genome_fasta "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz" \
  --transcriptome_gtf "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz" \
  --target "cellranger:bd_rhapsody" \
  --output_fasta "reference.fa.gz" \
  --output_gtf "reference.gtf.gz" \
  --output_cellranger "reference_cellranger.tar.gz" \
  --output_bd_rhapsody "reference_bd_rhapsody.tar.gz" \
  --subset_regex "chr1" \
  --publish_dir $OUT