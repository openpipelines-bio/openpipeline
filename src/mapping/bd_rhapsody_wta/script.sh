#!/bin/bash

# yq eval 'del(.. | select(has("dockerPull")))' cwl/v1.9.1/rhapsody_wta_1.9.1.cwl -j > rhapsody_wta_1.9.1_nodocker.cwl
# doesn't work

mkdir -p $par_output

cat > $par_output/config.yml << HERE
#!/usr/bin/env cwl-runner

cwl:tool: rhapsody

# This is a template YML file used to specify the inputs for a BD Genomics WTA Rhapsody Analysis pipeline run. See the
# BD Genomics Analysis Setup User Guide (Doc ID: 47383) for more details. Enter the following information:

## Reference_Genome (required) - Path to STAR index for tar.gz format. See Doc ID: 47383 for instructions to obtain pre-built STAR index file.
Reference_Genome:
   class: File
   location: "$(realpath $par_reference_genome)"

## Transcriptome_Annotation (required) - Path to GTF annotation file
Transcriptome_Annotation:
   class: File
   location: "$(realpath $par_transcriptome_annotation)"

## Reads (required) - Path to your read files in the FASTQ.GZ format. You may specify as many R1/R2 read pairs as you want.
Reads:
HERE

# process fastq files
IFS=:
set -f
for val in $par_input; do
  unset IFS
  cat >> $par_output/config.yml << HERE
 - class: File
   location: "$(realpath $val)"
HERE
done
set +f

# Add abseq reference, if specified
if [ ! -z "$par_abseq_reference" ]; then
  cat >> $par_output/config.yml << HERE

## AbSeq_Reference (optional) - Path to the AbSeq reference file in FASTA format.  Only needed if BD AbSeq Ab-Oligos are used.
AbSeq_Reference:
HERE

  # process abseq reference files
  IFS=:
  set -f
  for val in $par_abseq_reference; do
    unset IFS

    cat >> $par_output/config.yml << HERE
 - class: File
   location: "$(realpath $val)"
HERE
  done
  set +f
fi

# Add supplemental reference, if specified
if [ ! -z "$par_supplemental_reference" ]; then
  cat >> $par_output/config.yml << HERE

# Supplemental_Reference (optional) - Path to the supplemental reference file in FASTA format.  Only needed if there are additional transgene sequences used in the experiment.
Supplemental_Reference:
HERE

  # process supplemental reference files
  IFS=:
  set -f
  for val in $par_supplemental_reference; do
    unset IFS

    cat >> $par_output/config.yml << HERE
 - class: File
   location: "$(realpath $val)"
HERE
  done
  set +f
fi

# Add exact cell count, if specified
if [ ! -z "$par_exact_cell_count" ]; then
  cat >> $par_output/config.yml << HERE

## Exact cell count - Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count
Exact_Cell_Count: $par_exact_cell_count
HERE
fi

# Add Disable Refined Putative Cell Calling, if specified
if [ ! -z "$par_disable_putative_calling" ]; then
  cat >> $par_output/config.yml << HERE

## Disable Refined Putative Cell Calling - Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).
## The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.
## Does not apply if Exact Cell Count is set. Values can be true or false.
Basic_Algo_Only: $par_disable_putative_calling
HERE
fi

# add subsample, if specified
if [ ! -z "$par_subsample" ]; then
  cat >> $par_output/config.yml << HERE

## Subsample (optional) - A number >1 or fraction (0 < n < 1) to indicate the number or percentage of reads to subsample.
Subsample: $par_subsample
HERE
fi

cd $par_output

cwl-runner $resources_dir/rhapsody_wta_1.9.1_nodocker.cwl config.yml
