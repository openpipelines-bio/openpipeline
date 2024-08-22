#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=reference_gencodev41_chr1
OUT=resources_test/$ID

mkdir -p "$OUT" 

wget "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip" -O "$OUT/ERCC92.zip"

# Download JASPAR files for reference building
# Source of the code below: https://support.10xgenomics.com/single-cell-atac/software/release-notes/references#GRCh38-2020-A-2.0.0
motifs_url="https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"
motifs_in="${OUT}/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"

if [ ! -f "$motifs_in" ]; then
    curl -sS "$motifs_url" > "$motifs_in"
fi

# Change motif headers so the human-readable motif name precedes the motif
# identifier. So ">MA0004.1    Arnt" -> ">Arnt_MA0004.1".
motifs_modified="${OUT}/$(basename "$motifs_in").modified"
awk '{
    if ( substr($1, 1, 1) == ">" ) {
        print ">" $2 "_" substr($1,2)
    } else {
        print
    }
}' "$motifs_in" > "$motifs_modified"


cat > /tmp/params.yaml << HERE
param_list:
  - id: "$ID"
    genome_fasta: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz"
    transcriptome_gtf: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz"
    target: ["bd_rhapsody", "cellranger_arc"] 
    output_fasta: "reference.fa.gz"
    output_gtf: "reference.gtf.gz"
    non_nuclear_contigs: null
    output_cellranger_arc: "reference_cellranger.tar.gz"
    output_bd_rhapsody: "reference_bd_rhapsody.tar.gz"
    bdrhap_extra_star_params: "--genomeSAindexNbases 12 --genomeSAsparseD 2"
    motifs_file: "$motifs_modified"
    subset_regex: "chr1"
HERE

nextflow \
  run . \
  -main-script target/nextflow/workflows/ingestion/make_reference/main.nf \
  -profile docker \
  -c ./src/workflows/utils/labels_ci.config \
  -params-file /tmp/params.yaml \
  --publish_dir $OUT \
  -resume
