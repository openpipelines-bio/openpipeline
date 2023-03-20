#!/bin/bash

set -ex

echo ">>> Running executable"
$meta_executable \
    --bam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --fasta_reference "$meta_resources_dir/cellranger_tiny_fastq/cellranger_tiny_ref/fasta/genome.fa" \
    --output freebayes_result/ \
    --region chr21.part \
    --vcf mixed_variant.vcf

[[ ! -f freebayes_result/mixed_variant.vcf ]] && echo "Output VCF file could not be found!" && exit 1

echo ">>> Test finished successfully"
