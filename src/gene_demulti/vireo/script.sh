#!/bin/bash

set -eo pipefail

# add additional params
extra_params=( )

if [ ! -z "$par_vartrixData" ]; then 
  extra_params+=( "--vatrixData $par_vartrixData" )
fi
  
if [ ! -z "$par_donorFile" ]; then 
  extra_params+=( "--donorFile $par_donorFile" )
fi

if [ ! -z "$par_noDoublet" ]; then 
  extra_params+=( "--noDoublet" )
fi

if [ ! -z "$par_extraDonorMode" ] ; then 
  extra_params+=( "--extraDonorMode $par_extraDonorMode" )
fi

if [ ! -z "$par_forceLearnGT" ]; then 
  extra_params+=( "--forceLearnGT" )
fi

if [ ! -z "$par_ASEmode" ]; then 
  extra_params+=( "--ASEmode" )
fi

if [ ! -z "$par_noPlot" ]; then 
  extra_params+=( "--noPlot" )
fi

if [ ! -z "$par_randSeed" ]; then 
  extra_params+=( "--randSeed $par_randSeed" )
fi

if [ ! -z "$par_cellRange" ]; then 
  extra_params+=( "--cellRange $par_cellRange" )
fi

if [ ! -z "$par_callAmbientRNAs" ]; then 
  extra_params+=( "--callAmbientRNAs" )
fi

if [ ! -d "$par_output" ]; then
  mkdir $par_output
fi

vireo --cellData $par_cellData --nDonor $par_nDonor --genoTag $par_genoTag --nInit $par_nInit \
       --extraDonor $par_extraDonor --nproc $par_nproc --out ${par_output} ${extra_params[@]} 
