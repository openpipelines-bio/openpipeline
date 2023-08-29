#!/bin/bash

# so let's do it separately
rm images/concepts/fig_*.svg

for id in cell modality_rna modality_adt modality_vdj modality_atac workflow_multiomics_rna_singlesample workflow_multiomics_rna_multisample workflow_multiomics_adt_singlesample workflow_multiomics_adt_multisample; do
  inkscape --export-type="svg" --export-id="$id" --export-id-only images/concepts/fig.svg
  svgo images/concepts/fig_${id}.svg
done