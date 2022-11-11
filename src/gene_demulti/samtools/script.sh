#!/bin/bash
if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi

samtools view -S -b -q 10 -F 3844 $par_bam > ${par_output}filtered.bam
samtools index ${par_output}filtered.bam ${par_output}filtered.bam.bai
umi_tools dedup --stdin=filtered.bam --extract-umi-method=tag --umi-tag=UR --cell-tag=CB --log=logfile > ${par_output}no_dup.bam
samtools sort ${par_output}no_dup.bam -o ${par_output}sorted.bam
samtools index ${par_output}sorted.bam ${par_output}sorted.bam.bai
