#!/bin/bash

set -eou pipefail

## VIASH START
par_genome_fasta="resources_test/reference_gencodev41_chr1/reference.fa.gz"
par_transcriptome_gtf="resources_test/reference_gencodev41_chr1/reference.gtf.gz"
par_output="gencode_v41_annotation_cellranger.tar.gz"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# just to make sure
par_genome_fasta=`realpath $par_genome_fasta`
par_transcriptome_gtf=`realpath $par_transcriptome_gtf`
par_output=`realpath $par_output`


echo "> Unzipping input files"
unpigz -c "$par_genome_fasta" > "$tmpdir/genome.fa"

echo "> Building star index"
cd "$tmpdir"
cellranger mkref \
  --fasta "$tmpdir/genome.fa" \
  --genes "$par_transcriptome_gtf" \
  --genome output \
  ${par_reference_version:+--ref-version $par_reference_version} \
  ${meta_cpus:+--nthreads $meta_cpus} \
  ${meta_memory_gb:+--memgb $(($meta_memory_gb-2))} # always keep 2 gb for the OS itseld

echo "> Creating archive"
tar --use-compress-program="pigz -k " -cf "$par_output" -C "$tmpdir/output" .