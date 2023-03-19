#!/bin/bash

set -eo pipefail

# add additional params
extra_params=( )
  
if [ ! -z "$par_restarts" ]; then
  extra_params+=( "--restarts $par_restarts" )
fi

if [ ! -z "$par_common_variants" ]; then
  extra_params+=( "--common_variants $par_common_variants" )
fi

if [ ! -z "$par_known_genotypes" ] ; then
  extra_params+=( "--known_genotypes $par_known_genotypes" )
fi

if [ ! -z "$par_known_genotypes_sample_names" ] ; then
  extra_params+=( "--known_genotypes_sample_names $par_known_genotypes_sample_names" )
fi

if [ "$par_skip_remap" = true ]; then
  extra_params+=( "--skip_remap True" )
fi

if [ "$par_ignore" = true ]; then
  extra_params+=( "--ignore True" )
fi

if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi

/opt/souporcell/souporcell_pipeline.py --bam $par_bam --fasta $par_fasta --barcodes $par_barcodes \
--clusters $par_clusters --ploidy $par_ploidy --min_alt $par_min_alt --min_ref $par_min_ref \
--max_loci $par_max_loci --threads $par_threads --out_dir $par_output ${extra_params[@]}
cut -d$'\t' -f 1-2 ${par_output}clusters.tsv > ${par_output}assignment.tsv
sed 1d ${par_output}/assignment.tsv > ${par_output}res.tsv
{ echo cell$'\t'donor_id; cat ${par_output}res.tsv; } > ${par_output}assignment.tsv
rm ${par_output}res.tsv
