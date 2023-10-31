#!/bin/bash

set -eo pipefail

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

scSplit count \
  --vcf $par_vcf \
  --bam $par_bam \
  --bar $par_bar \
  --tag $par_tag \
  --ref $par_ref \
  --alt $par_alt \
  --out $par_output \
  ${par_com:+--com $par_com}

scSplit run \
  --ref "$par_output/$par_ref" \
  --alt "$par_output/$par_alt" \
  --out $par_output \
  --num $par_num \
  ${par_sub:+--sub $par_sub} \
  ${par_ems:+--ems $par_ems} \
  ${par_dbl:+--dbl $par_dbl} \
  ${par_vcf_known:+--vcf $par_vcf_known}

if [ "$par_geno" = true ]; then
  scSplit genotype \
    --ref "$par_output/$par_ref" \
    --alt "$par_output/$par_alt" \
    --psc "$par_output/$par_psc" \
    "$par_output"
fi

echo "cell,donor_id" > "$par_output/cell_annotation.csv"
sed -e '1d' -e 's/SNG-//g' "$par_output/scSplit_result.csv" | 
sed 's/\t/,/g' | awk 'BEGIN{FS=OFS=","} { if ($2 ~ /^DBL-/) $2 = "doublet"; print }' \
>> "$par_output/cell_annotation.csv"
