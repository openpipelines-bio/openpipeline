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

nextflow \
  run . \
  -main-script target/nextflow/workflows/ingestion/make_reference/main.nf \
  -profile docker \
  --id "$ID" \
  -c ./src/workflows/utils/labels_ci.config \
  --genome_fasta "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz" \
  --transcriptome_gtf "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz" \
  --target "cellranger;bd_rhapsody" \
  --output_fasta "reference.fa.gz" \
  --output_gtf "reference.gtf.gz" \
  --output_cellranger "reference_cellranger.tar.gz" \
  --output_bd_rhapsody "reference_bd_rhapsody.tar.gz" \
  --subset_regex "chr1" \
  --publish_dir $OUT \
  -resume


viash run src/reference/build_bdrhap2_reference/config.vsh.yaml -- \
  --genome_fasta "$OUT/reference.fa.gz" \
  --gtf "$OUT/reference.gtf.gz" \
  --reference_archive "$OUT/reference_bd_rhapsody_v2.tar.gz" \
  --extra_star_params "--genomeSAindexNbases 4" \
  ---cpus 2 

# Generate data processed by BD Rhapsody 2.x.x
bdrhap_5kjrt_resources="resources_test/bdrhap_5kjrt"

viash run src/mapping/bd_rhapsody2/config.vsh.yaml -- \
  --reads "$bdrhap_5kjrt_resources/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  --reads "$bdrhap_5kjrt_resources/raw/12WTA_S1_L432_R2_001_subset.fastq.gz" \
  --reads "$bdrhap_5kjrt_resources/raw/12ABC_S1_L432_R1_001_subset.fastq.gz" \
  --reads "$bdrhap_5kjrt_resources/raw/12ABC_S1_L432_R2_001_subset.fastq.gz" \
  --reference_archive "$OUT/reference_bd_rhapsody_v2.tar.gz" \
  --abseq_reference "$bdrhap_5kjrt_resources/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta" \
  --output_dir "$bdrhap_5kjrt_resources/processed2" \
  --cell_calling_data "mRNA" \
  --exact_cell_count 4900 \
  ---memory 10gb \
  ---cpus 2
