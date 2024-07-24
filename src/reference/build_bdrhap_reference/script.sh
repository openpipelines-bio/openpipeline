#!/bin/bash

set -eou pipefail

## VIASH START
par_genome_fasta="temp/reference.fa.gz"
par_transcriptome_gtf="temp/reference.gtf.gz"
par_output="temp/reference_star.tar.gz"
meta_cpus=20
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

meta_cpus="${meta_cpus:-1}"

# process params
extra_params=( )

if [ ! -z "$meta_cpus" ]; then 
  extra_params+=( "--runThreadN $meta_cpus" )
fi

echo "> Unzipping input files"
unpigz -c "$par_genome_fasta" > "$tmpdir/genome.fa"
unpigz -c "$par_transcriptome_gtf" > "$tmpdir/transcriptome.gtf"

echo "> Building star index"
mkdir "$tmpdir/out"
STAR \
  --runMode genomeGenerate \
  --genomeDir "$tmpdir/out" \
  --genomeFastaFiles "$tmpdir/genome.fa" \
  --sjdbGTFfile "$tmpdir/transcriptome.gtf" \
  --sjdbOverhang 100 \
  --genomeSAindexNbases 11 \
  "${extra_params[@]}"

echo "> Creating archive"
tar --use-compress-program="pigz -k " -cf "$par_output" -C "$tmpdir/out" .