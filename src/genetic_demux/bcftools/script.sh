#!/bin/bash
if [ ! -d "$par_output" ]; then
  mkdir -p $par_output
fi

IFS=";" read -a vcf_list <<< $par_vcf


if [ "$par_concat" = true ] && [ "$par_filter" = true ] ; then
  bcftools concat -o "$par_output/concated_chroms.vcf" ${vcf_list[@]}
  bcftools sort "$par_output/concated_chroms.vcf" -o "$par_output/sorted_concated_chroms.vcf"
  bcftools filter -i "QUAL>$par_filter_qual" "$par_output/sorted_concated_chroms.vcf" -o "$par_output/filtered_sorted_concated_chroms.vcf"
  
elif [ "$par_filter" = true ] ; then
  bcftools filter -i "QUAL>$par_filter_qual" ${vcf_list[@]} -o "$par_output/filtered.vcf"
    
else
  bcftools concat -o "$par_output/concated_chroms.vcf" ${vcf_list[@]}
  bcftools sort "$par_output/concated_chroms.vcf" -o "$par_output/sorted_concated_chroms.vcf"
fi
