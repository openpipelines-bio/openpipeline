# set -eo pipefail

## VIASH START
meta_resources_dir="./resources_test"
meta_executable="./target/docker/mapping/star_build_reference/star_build_reference"
## VIASH END

"$meta_executable" \
  --genome_fasta "$meta_resources_dir/cellranger_tiny_fastq/cellranger_tiny_ref/fasta/genome.fa" \
  --transcriptome_gtf "$meta_resources_dir/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz" \
  --output "./star_reference_test" \
  --genomeSAindexNbases 7 \
  ---cpus 8

if [ ! -f "./star_reference_test/Genome" ]; then
    echo "Genome file could not be found in the output directory";
    exit 1
fi