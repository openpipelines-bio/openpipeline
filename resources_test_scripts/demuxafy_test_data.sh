#!/bin/bash

set -eo pipefail


# settings
ID=demuxafy_test_data
OUT=resources_test/$ID
DIR="$OUT"

mkdir -p "$OUT"
cd "$OUT"
# download demuxafy test dataset
wget https://www.dropbox.com/s/m8u61jn4i1mcktp/TestData4PipelineSmall.tar.gz
tar -xf TestData4PipelineSmall.tar.gz
# bam and vcf file
cp TestData4PipelineSmall/test_dataset.vcf .
# extract chr from vcf file
grep -w '^#\|^#CHROM\|^[1-2]' test_dataset.vcf > test_dataset_chr1_2.vcf
grep -w '^#\|^#CHROM\|^[3-4]' test_dataset.vcf > test_dataset_chr3_4.vcf

# barcode list
cp TestData4PipelineSmall/test_dataset/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/barcodes.tsv .

# subsetted bam and bai for souporcell
wget https://www.dropbox.com/s/7ew5lt0msf4z5gj/chr_1_pooled.sorted.bam
wget https://www.dropbox.com/s/tpplbj9sab9b2p4/chr_1_pooled.sorted.bam.bai

# variants from mixed sample
wget https://www.dropbox.com/s/btir7ge4kzc7tu1/mixed_variant.vcf

# dsc_pileup output
wget https://www.dropbox.com/s/17hj9i0yavtezx1/dsc_pileup.zip
unzip dsc_pileup.zip

# subsetted human genome reference
wget https://www.dropbox.com/s/ynlce3g7nwxthwg/genome_chr1.fa

# remove unnecessary files
rm -rf TestData4PipelineSmall
rm TestData4PipelineSmall.tar.gz
rm dsc_pileup.zip