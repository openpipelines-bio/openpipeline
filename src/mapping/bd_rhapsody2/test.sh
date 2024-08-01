#!/bin/bash

set -e

## VIASH START
meta_executable="target/docker/mapping/bd_rhapsody2/bd_rhapsody_2"
meta_resources_dir="src/mapping/bd_rhapsody2"
## VIASH END

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_doesnt_exist() {
  [ ! -f "$1" ] || { echo "File '$1' exists but shouldn't" && exit 1; }
}
assert_file_empty() {
  [ ! -s "$1" ] || { echo "File '$1' is not empty but should be" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains() {
  grep -q "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}
assert_file_contains_regex() {
  grep -q -E "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_file_not_contains_regex() {
  grep -q -E "$2" "$1" && { echo "File '$1' contains '$2' but shouldn't" && exit 1; }
}

#########################################################################################

# generate index. this is not part of the unit test, but data we need to run
# the bd rhapsody pipeline.
# cwl_file="$meta_resources_dir/make_rhap_reference_2.2.1_nodocker.cwl"
cwl_file="src/mapping/bd_rhapsody2/make_rhap_reference_2.2.1_nodocker.cwl"
reference_file="resources_test/reference_gencodev41_chr1/Rhap_reference.tar.gz"
# genome_fasta="$meta_resources_dir/reference.fa.gz"
# genome_gtf="$meta_resources_dir/reference.gtf.gz"
genome_fasta="resources_test/reference_gencodev41_chr1/reference.fa.gz"
genome_gtf="resources_test/reference_gencodev41_chr1/reference.gtf.gz"

# echo "============="
# gunzip "$genome_fasta" > "$meta_resources_dir/reference.fa" 
# gunzip "$genome_gtf" > "$meta_resources_dir/reference.gtf" 
# echo "++++++++++++++"

cwl-runner \
  --no-container \
  --preserve-entire-environment \
  --outdir "resources_test/reference_gencodev41_chr1" \
  "$cwl_file" \
  --Genome_fasta "$genome_fasta" \
  --Gtf "$genome_gtf" \
  --Extra_STAR_params "--genomeSAindexNbases 4"

#############################################

# # viash run src/mapping/bd_rhapsody2/config.vsh.yaml -- \
# #   --reads "resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
# #   --reads "resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz" \
# #   --reference_archive "../biobox/reference_gencodev41_chr1.tar.gz" \
# #   --output_dir output_large \
# #   --cell_calling_data mRNA \
# #   --exact_cell_count 4900 \
# #   ---memory 15gb \
# #   ---cpus 10

# echo ">> Run $meta_name"
# "$meta_executable" \
#   --reads "$meta_resources_dir/12WTA_S1_L432_R1_001_subset.fastq.gz" \
#   --reads "$meta_resources_dir/12WTA_S1_L432_R2_001_subset.fastq.gz" \
#   --reference_archive "$REFERENCE_FILE" \
#   --output_dir output \
#   ${meta_cpus:+---cpus $meta_cpus} \
#   ${meta_memory_mb:+---memory ${meta_memory_mb}MB} \
#   --reads ABCreads_R1.fq.gz \
#   --reads ABCreads_R2.fq.gz \
#   --abseq_reference bdabseq_smallpanel.fasta \
#   --exact_cell_count 2 \
#   --exclude_intronic_reads false \

# echo ">> Check if output exists"
# assert_file_exists "output/sample_Bioproduct_Stats.csv"
# assert_file_exists "output/sample_RSEC_MolsPerCell_Unfiltered_MEX.zip"
# assert_file_exists "output/Logs"
# assert_file_exists "output/sample_Metrics_Summary.csv"

# # echo ">> Check if output contents are not empty"
# assert_file_not_empty "output/sample_Bioproduct_Stats.csv"
# assert_file_not_empty "output/sample_RSEC_MolsPerCell_Unfiltered_MEX.zip"
# assert_file_not_empty "output/Logs"
# assert_file_not_empty "output/sample_Metrics_Summary.csv"

# # echo ">> Check if output contents are correct"
# # assert_file_contains "log.txt" "Number of input reads \\|	2"
# # assert_file_contains "log.txt" "Number of reads unmapped: too short \\|	1"
# # assert_file_contains "log.txt" "Uniquely mapped reads number \\|	1"

# #########################################################################################

# # TODO: add test with ABC, VDJ, SMK, and ATAC

# echo "> Test successful"
