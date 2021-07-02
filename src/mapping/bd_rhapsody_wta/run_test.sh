#!/bin/bash

mkdir output

./bd_rhapsody_wta \
  -i sample_R1_.fastq.gz \
  -i sample_R2_.fastq.gz \
  -r GRCh38-PhiX-gencodev29-20181205.tar.gz \
  -t gencodev29-20181205.gtf \
  --subsample 0.5 \
  -o output/

[[ ! -f output/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1


: << 'COMMENT'

mkdir reference
cd reference/

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz

for chr in {21,22}; do
  mkdir chr$chr  
  awk '/>chr/{p=0} />chr'$chr' /{p=1} p' GRCh38.primary_assembly.genome.fa > chr$chr/GRCh38.primary_assembly.genome.chr$chr.fa
  grep chr$chr gencode.v29.annotation.gtf > chr$chr/gencode.v29.annotation.chr$chr.gtf  
  STAR --runThreadN 20 --runMode genomeGenerate --genomeDir STAR_chr$chr --genomeFastaFiles chr$chr/GRCh38.primary_assembly.genome.chr$chr.fa --sjdbGTFfile chr$chr/gencode.v29.annotation.chr$chr.gtf --sjdbOverhang 100
  tar -czvf chr$chr/GRCh38.primary_assembly.genome.index.chr$chr.tar.gz STAR_chr$chr
  rm -r STAR_chr$chr
done

./bd_rhapsody_wta \
  -i sample_R1_.fastq.gz \
  -i sample_R2_.fastq.gz \
  -r reference/chr21/GRCh38.primary_assembly.genome.index.chr21.tar.gz \
  -t reference/chr21/gencode.v29.annotation.chr21.gtf \
  --subsample 0.5 \
  -o output/
  
./bd_rhapsody_wta \
  -i sample_R1_.fastq.gz \
  -i sample_R2_.fastq.gz \
  -r reference/chr21/STAR_chr21.tar.gz \
  -t reference/chr21/gencode.v29.annotation.chr21.gtf \
  --subsample 0.5 \
  -o output/
COMMENT

echo ">>> Test finished successfully"
