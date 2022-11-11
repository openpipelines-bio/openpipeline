#!/bin/bash

set -eo pipefail

# add additional params
extra_params=( )

if [ ! -z "$par_samFile" ]; then
  extra_params+=( "--samFile $par_samFile" )
fi
  
if [ ! -z "$par_samFileList" ]; then
  extra_params+=( "--samFileList $par_samFileList" )
fi

if [ ! -z "$par_regionsVCF" ]; then
  extra_params+=( "--regionsVCF $par_regionsVCF" )
fi

if [ ! -z "$par_targetsVCF" ] ; then
  extra_params+=( "--targetsVCF $par_targetsVCF" )
fi

if [ ! -z "$par_barcodeFile" ]; then
  extra_params+=( "--barcodeFile $par_barcodeFile" )
fi

if [ ! -z "$par_sampleList" ]; then
  extra_params+=( "--sampleList $par_sampleList" )
fi

if [ ! -z "$par_sampleIDs" ]; then
  extra_params+=( "--sampleIDs $par_sampleIDs" )
fi

if [ "$par_genotype" = true ]; then
  extra_params+=( "--genotype" )
fi

if [ "$par_gzip" = true ]; then
  extra_params+=( "--gzip" )
fi

if [ "$par_printSkipSNPs" = true ]; then
  extra_params+=( "--printSkipSNPs" )
fi

if [ ! -z "$par_chrom" ]; then
  extra_params+=( "--chrom $par_chrom" )
fi

if [ "$par_doubletGL" = true ]; then
  extra_params+=( "--doubletGL")
fi

if [ ! -z "$par_inclFLAG" ]; then
  extra_params+=( "--inclFLAG $par_inclFLAG" )
fi

if [ ! -z "$par_exclFLAG" ]; then
  extra_params+=( "--exclFLAG $par_exclFLAG" )
fi

if [ "$par_countORPHAN" = true ]; then
  extra_params+=( "--countORPHAN" )
fi
          
cellsnp-lite --nproc $par_nproc --cellTAG $par_cellTAG --UMItag $par_UMItag \
             --minCOUNT $par_minCOUNT --minMAF $par_minMAF --minLEN $par_minLEN \
             --minMAPQ $par_minMAPQ --outDir $par_output ${extra_params[@]}