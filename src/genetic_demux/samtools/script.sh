#!/bin/bash
if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi
samtools view -S -b -q 10 -F 3844 "$par_bam" > "${par_output}/filtered.bam"
cd $par_output
samtools index filtered.bam filtered.bam.bai
umi_tools dedup --stdin=filtered.bam --extract-umi-method=tag --umi-tag=UR --cell-tag=CB --log=logfile > no_dup.bam
samtools sort no_dup.bam -o sorted.bam
samtools index sorted.bam sorted.bam.bai
