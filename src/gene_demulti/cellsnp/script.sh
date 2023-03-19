#!/bin/bash

set -eo pipefail

# Unset flags if they equal 'false'
[[ "$par_genotype" == "false" ]] && unset par_genotype
[[ "$par_gzip" == "false" ]] && unset par_gzip
[[ "$par_printSkipSNPs" == "false" ]] && unset par_printSkipSNPs
[[ "$par_doubletGL" == "false" ]] && unset par_doubletGL
[[ "$par_countORPHAN" == "false" ]] && unset par_countORPHAN

cellsnp-lite \
  ${meta_cpus:+--nproc $meta_cpus} \
  --cellTAG $par_cellTAG \
  --UMItag $par_UMItag \
  --minCOUNT $par_minCOUNT \
  --minMAF $par_minMAF \
  --minLEN $par_minLEN \
  --minMAPQ $par_minMAPQ \
  --outDir $par_output \
  ${par_samFile:+--samFile $par_samFile} \
  ${par_samFileList:+--samFileList $par_samFileList} \
  ${par_regionsVCF:+--regionsVCF $par_regionsVCF} \
  ${par_targetsVCF:+--targetsVCF $par_targetsVCF} \
  ${par_barcodeFile:+--barcodeFile $par_barcodeFile} \
  ${par_sampleList:+--sampleList $par_sampleList} \
  ${par_sampleIDs:+--sampleIDs $par_sampleIDs} \
  ${par_genotype:+--genotype} \
  ${par_gzip:+--gzip} \
  ${par_printSkipSNPs:+--printSkipSNPs} \
  ${par_chrom:+--chrom $par_chrom} \
  ${par_doubletGL:+--doubletGL} \
  ${par_inclFLAG:+--inclFLAG $par_inclFLAG} \
  ${par_exclFLAG:+--exclFLAG $par_exclFLAG} \
  ${par_countORPHAN:+--countORPHAN}
