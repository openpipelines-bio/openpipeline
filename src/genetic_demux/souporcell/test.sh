#!/bin/bash

set -exo pipefail

data="$meta_resources_dir/demuxafy_test_data"

echo ">>> Running executable"
$meta_executable \
    ---cpus "${meta_cpus:-1}" \
    --bam "$data/chr_1_pooled.sorted.bam" \
    --bam_index "$data/chr_1_pooled.sorted.bam.bai" \
    --fasta "$data/genome_chr1.fa" \
    --barcodes "$data/barcodes.tsv" \
    --common_variants "$data/test_dataset.vcf" \
    --clusters 14 \
    --skip_remap \
    --output soup_result/

[[ ! -f soup_result/clusters.tsv ]] && echo "Output donor assignment file could not be found!" && exit 1
[[ ! -f soup_result/cell_annotation.csv ]] && echo "Output cell annotation file could not be found!" && exit 1

echo ">>> Running executable with khm clustering and souporcell3"
$meta_executable \
    ---cpus "${meta_cpus:-1}" \
    --bam "$data/chr_1_pooled.sorted.bam" \
    --bam_index "$data/chr_1_pooled.sorted.bam.bai" \
    --fasta "$data/genome_chr1.fa" \
    --barcodes "$data/barcodes.tsv" \
    --common_variants "$data/test_dataset.vcf" \
    --clusters 14 \
    --clustering_method khm \
    --souporcell3 \
    --restarts 10 \
    --skip_remap \
    --output soup_result_khm/ 2>&1 | tee soup_result_khm.log

[[ ! -f soup_result_khm/cell_annotation.csv ]] && echo "Output cell annotation file could not be found!" && exit 1
grep -q -- "--restarts 10 .* -s true -m khm" soup_result_khm.log || { echo "Clustering options were not passed to souporcell!"; exit 1; }

echo ">>> Running executable with known genotypes"
sample_names=$(grep -m1 '^#CHROM' "$data/test_dataset.vcf" | cut -f10- | tr '\t' ';')
$meta_executable \
    ---cpus "${meta_cpus:-1}" \
    --bam "$data/chr_1_pooled.sorted.bam" \
    --bam_index "$data/chr_1_pooled.sorted.bam.bai" \
    --fasta "$data/genome_chr1.fa" \
    --barcodes "$data/barcodes.tsv" \
    --known_genotypes "$data/test_dataset.vcf" \
    --known_genotypes_sample_names "$sample_names" \
    --clusters 14 \
    --skip_remap \
    --output soup_result_known/ 2>&1 | tee soup_result_known.log

[[ ! -f soup_result_known/cell_annotation.csv ]] && echo "Output cell annotation file could not be found!" && exit 1
grep -q -- "--known_genotypes_sample_names ${sample_names//;/ }" soup_result_known.log || { echo "Sample names were not passed to souporcell!"; exit 1; }

echo ">>> Test finished successfully"
