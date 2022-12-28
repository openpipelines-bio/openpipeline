#!/bin/bash

set -eo pipefail

# add additional params
extra_params=( )
extra_params_count=( )
if [ ! -z "$par_sub" ]; then
  extra_params+=( "--sub $par_sub" )
fi
  
if [ ! -z "$par_ems" ]; then
  extra_params+=( "--ems $par_ems" )
fi

if [ ! -z "$par_dbl" ]; then
  extra_params+=( "--dbl $par_dbl" )
fi

if [ ! -z "$par_vcf_known" ] ; then
  extra_params+=( "--vcf $par_vcf_known" )
fi

if [ ! -z "$par_com" ] ; then
  extra_params_count+=( "--com $par_com" )
fi

if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi

scSplit count --vcf $par_vcf --bam $par_bam \
        --bar $par_bar --tag $par_tag --ref $par_ref --alt $par_alt \
        --out $par_output ${extra_params_count[@]}
scSplit run --ref ${par_output}$par_ref \
        --alt ${par_output}$par_alt --out $par_output --num $par_num ${extra_params[@]}

if [ "$par_geno" = true ]; then
    scSplit genotype --ref ${par_output}$par_ref \
        --alt ${par_output}$par_alt --psc ${par_output}${par_psc} $par_output
fi

