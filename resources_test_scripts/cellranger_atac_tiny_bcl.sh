#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# settings
ID=cellranger_atac_tiny_bcl
OUT="resources_test/$ID/"
DIR="$OUT"
REFERENCE_DIR=resources_test/reference_gencodev41_chr1

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

viash ns build -q "download_file|cellranger_atac_mkfastq" -p docker --setup cb

# download bcl data
if [ ! -f "${OUT}/bcl/sample_sheet.csv" ]; then
  mkdir -p "$OUT/bcl"

  # download tar gz
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/supp/cell-atac/cellranger-atac-tiny-bcl-1.0.0.tar.gz \
    --output "${OUT}/bcl/cellranger-atac-tiny-bcl-1.0.0.tar.gz"
  
  # untar
  
  tar -xf "${OUT}/bcl/cellranger-atac-tiny-bcl-1.0.0.tar.gz" \
    --strip-components=1 \
    -C "$OUT/bcl"

  # remove tar
  rm "${OUT}/bcl/cellranger-atac-tiny-bcl-1.0.0.tar.gz"

  # Download the layout file. It contains info about the samples (1 in this case) and lanes
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/supp/cell-atac/cellranger-atac-tiny-bcl-simple-1.0.0.csv \
    --output "${OUT}/bcl/layout.csv"

  # download sample sheet
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/supp/cell-atac/cellranger-atac-tiny-bcl-samplesheet-1.0.0.csv \
    --output "${OUT}/bcl/sample_sheet.csv"
fi

# Download JASPAR files for reference building
# Source of the code below: https://support.10xgenomics.com/single-cell-atac/software/release-notes/references#GRCh38-2020-A-2.0.0
motifs_url="https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"
motifs_in="${REFERENCE_DIR}/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"

if [ ! -f "$motifs_in" ]; then
    curl -sS "$motifs_url" > "$motifs_in"
fi

# Change motif headers so the human-readable motif name precedes the motif
# identifier. So ">MA0004.1    Arnt" -> ">Arnt_MA0004.1".
motifs_modified="${REFERENCE_DIR}/$(basename "$motifs_in").modified"
awk '{
    if ( substr($1, 1, 1) == ">" ) {
        print ">" $2 "_" substr($1,2)
    } else {
        print
    }
}' "$motifs_in" > "$motifs_modified"

if [ ! -d "${OUT}/fastqs" ]; then
  mkdir -p "$OUT/fastqs"

  target/docker/demux/cellranger_atac_mkfastq/cellranger_atac_mkfastq \
    --input "${OUT}/bcl" \
    --csv "${OUT}/bcl/layout.csv" \
    --output "${OUT}/fastqs"
fi

# Create reference
target/docker/reference/build_cellranger_arc_reference \
  --genome_fasta "${REFERENCE_DIR}/reference.fa.gz" \
  --annotation_gtf "${REFERENCE_DIR}/reference.gtf.gz" \
  --motifs_file "${REFERENCE_DIR}/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt.modified" \
  --output "${REFERENCE_DIR}/reference_cellranger_arc.tar.gz" \
  --organism "Homo_sapiens" \
  --genome "GRCh38"

# Create count matrices
if [ ! -d "${OUT}/counts" ]; then
  mkdir -p "$OUT/counts"
  
  target/docker/mapping/cellranger_atac_count \
    --input "${OUT}/fastqs/HJN3KBCX2/test_sample/" \
    --reference "${REFERENCE_DIR}/reference_cellranger_arc.tar.gz" \
    --output "${OUT}/counts"
