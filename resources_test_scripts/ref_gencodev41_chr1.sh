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
  -main-script target/nextflow/workflows/ingestion/make_reference/main.nf \
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


cwl_file="src/mapping/bd_rhapsody2/make_rhap_reference_2.2.1_nodocker.cwl"
reference_small_gtf="$OUT/reference.gtf"
reference_small_fa="$OUT/reference.fa"

viash run src/reference/build_bdrhap2_reference/config.vsh.yaml -- \
  --genome_fasta "resources_test/reference_gencodev41_chr1/reference.fa.gz" \
  --gtf "resources_test/reference_gencodev41_chr1/reference.gtf.gz" \
  --reference_archive "Rhap_reference.tar.gz" \
  --extra_star_params "--genomeSAindexNbases 4" \
  ---cpus 2

mv Rhap_reference.tar.gz $OUT

viash run src/mapping/bd_rhapsody2/config.vsh.yaml -- \       
  --reads "resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  --reads "resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz" \
  --reference_archive "resources_test/reference_gencodev41_chr1/Rhap_reference.tar.gz" \
  --output_dir "processed2" \
  --cell_calling_data "mRNA" \
  --exact_cell_count 4900

mv processed2 resources_test/bdrhap_5kjrt
