#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"


export NXF_VER=21.04.1

PAR_GENOME="/scratch/GRCh38-PhiX-gencodev29-20181205.tar.gz"
PAR_TRANSCRIPTOME="/scratch/gencodev29-20181205.gtf"
PAR_INPUT='/app/project/CS000182/210729_BD_Rhapsody_NOPRODHPB0013_FNA/FASTQ/GC113427_AGCGTAGC-TATAGCCT.*.R[12].fastq.gz'

nextflow \
  run . \
  -main-script workflows/1_ingestion/bd_rhapsody_wta/main.nf \
  -entry auto_wf \
  --input "$PAR_INPUT" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  -resume \
  -w /scratch/nxf_work_openpipeline
