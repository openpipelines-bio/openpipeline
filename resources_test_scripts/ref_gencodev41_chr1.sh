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


# create small reference by subsetting to small region
START=30000
END=31500
CHR=chr1

seqkit grep -r -p "^$CHR\$" "$OUT/reference.fa.gz" | \
  seqkit subseq -r "$START:$END" | gzip > $OUT/reference_small.fa.gz

zcat "$OUT/reference.gtf.gz" | \
  awk -v FS='\t' -v OFS='\t' "
    \$1 == \"$CHR\" && \$4 >= $START && \$5 <= $END {
      \$4 = \$4 - $START + 1;
      \$5 = \$5 - $START + 1;
      print;
    }" | gzip > $OUT/reference_small.gtf.gz

# Build references for BD Rhapsody 2.x.x
viash run src/reference/build_bdrhap2_reference/config.vsh.yaml -- \
  --genome_fasta "$OUT/reference.fa.gz" \
  --gtf "$OUT/reference.gtf.gz" \
  --reference_archive "$OUT/Rhap_reference_full.tar.gz" \
  --extra_star_params "--genomeSAindexNbases 4" \
  ---cpus 2 
  
  viash run src/reference/build_bdrhap2_reference/config.vsh.yaml -- \
  --genome_fasta "$OUT/reference_small.fa.gz" \
  --gtf "$OUT/reference_small.gtf.gz" \
  --reference_archive "$OUT/Rhap_reference_small.tar.gz" \
  --extra_star_params "--genomeSAindexNbases 4" \
  ---cpus 2 

# Generate data processed by BD Rhapsody 2.x.x
bdrhap_5kjrt_resources="resources_test/bdrhap_5kjrt"

viash run src/mapping/bd_rhapsody2/config.vsh.yaml -- \
  --reads "$bdrhap_5kjrt_resources/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  --reads "$bdrhap_5kjrt_resources/raw/12WTA_S1_L432_R2_001_subset.fastq.gz" \
  --reads "$bdrhap_5kjrt_resources/raw/12ABC_S1_L432_R1_001_subset.fastq.gz" \
  --reads "$bdrhap_5kjrt_resources/raw/12ABC_S1_L432_R2_001_subset.fastq.gz" \
  --reference_archive "$OUT/Rhap_reference_small.tar.gz" \
  --abseq_reference "$bdrhap_5kjrt_resources/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta" \
  --output_dir "$bdrhap_5kjrt_resources/processed_small" \
  --cell_calling_data "mRNA" \
  --exact_cell_count 4900 \
  ---memory 10gb \
  ---cpus 2

rm "$OUT/Rhap_reference_full.tar.gz"