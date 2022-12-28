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
cp TestData4PipelineSmall/test_dataset/outs/pooled.sorted.bam.bai .
cp TestData4PipelineSmall/test_dataset/outs/pooled.sorted.bam .
cp TestData4PipelineSmall/test_dataset.vcf .
# extract chr from vcf file
grep -w '^#\|^#CHROM\|^[1-2]' test_dataset.vcf > test_dataset_chr1_2.vcf
grep -w '^#\|^#CHROM\|^[3-4]' test_dataset.vcf > test_dataset_chr3_4.vcf

# barcode list
cp TestData4PipelineSmall/test_dataset/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/barcodes.tsv .

# variants from mixed sample
wget https://www.dropbox.com/s/btir7ge4kzc7tu1/mixed_variant.vcf

# human genome reference
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
tar -xzvf refdata-cellranger-GRCh38-3.0.0.tar.gz
mv refdata-cellranger-GRCh38-3.0.0/fasta/* .

# remove unnecessary files
rm -rf TestData4PipelineSmall
rm -rf refdata-cellranger-GRCh38-3.0.0
rm TestData4PipelineSmall.tar.gz
rm refdata-cellranger-GRCh38-3.0.0.tar.gz
