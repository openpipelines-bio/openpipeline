#!/bin/bash

set -eo pipefail

# Unset flags if they equal 'false'
[[ "$par_genotype" == "false" ]] && unset par_genotype
[[ "$par_gzip" == "false" ]] && unset par_gzip
[[ "$par_print_skip_snps" == "false" ]] && unset par_print_skip_snps
[[ "$par_doublet_gl" == "false" ]] && unset par_doublet_gl
[[ "$par_count_orphan" == "false" ]] && unset par_count_orphan

cellsnp-lite \
  ${meta_cpus:+--nproc $meta_cpus} \
  --cellTAG $par_cell_tag \
  --UMItag $par_umi_tag \
  --minCOUNT $par_min_count \
  --minMAF $par_min_maf \
  --minLEN $par_min_len \
  --minMAPQ $par_min_mapq \
  --outDir $par_output \
  ${par_sam_file:+--samFile $par_sam_file} \
  ${par_sam_fileList:+--samFileList $par_sam_fileList} \
  ${par_regions_vcf:+--regionsVCF $par_regions_vcf} \
  ${par_targets_vcf:+--targetsVCF $par_targets_vcf} \
  ${par_barcode_file:+--barcodeFile $par_barcode_file} \
  ${par_sample_list:+--sampleList $par_sample_list} \
  ${par_sample_ids:+--sampleIDs $par_sample_ids} \
  ${par_genotype:+--genotype} \
  ${par_gzip:+--gzip} \
  ${par_print_skip_snps:+--printSkipSNPs} \
  ${par_chrom:+--chrom $par_chrom} \
  ${par_doublet_gl:+--doubletGL} \
  ${par_incl_flag:+--inclFLAG $par_incl_flag} \
  ${par_excl_flag:+--exclFLAG $par_excl_flag} \
  ${par_count_orphan:+--countORPHAN}
