#!/bin/bash

set -eo pipefail

## VIASH START
par_input="resources_test/cellranger_tiny_fastq/bam/possorted_genome_bam.bam"
par_transcriptome="resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf"
par_output="./foo/output.loom"
par_barcode=""
## VIASH END

extra_params=( )

if [ ! -z "$par_barcode" ]; then 
  extra_params+=( "--bcfile=$par_barcode" )
fi
if [ "$par_without_umi" == "true" ]; then
  extra_params+=( "--without-umi" )
fi
if [ ! -z "$meta_cpus" ]; then
  extra_params+=( "--samtools-threads" "$meta_cpus" )
fi
if [ ! -z "$meta_memory_mb" ]; then
  extra_params+=( "--samtools-memory" "$meta_memory_mb" )
fi

output_dir=`dirname "$par_output"`
sample_id=`basename "$par_output" .loom`

if (file `readlink -f "$par_transcriptome"` | grep -q compressed ) ; then
  # create temporary directory
  tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
  function clean_up {
      rm -rf "$tmpdir"
  }
  trap clean_up EXIT

  zcat "$par_transcriptome" > "$tmpdir/genes.gtf"
  par_transcriptome="$tmpdir/genes.gtf"
fi

velocyto run \
  "$par_input" \
  "$par_transcriptome" \
  "${extra_params[@]}" \
  --outputfolder "$output_dir" \
  --sampleid "$sample_id"
