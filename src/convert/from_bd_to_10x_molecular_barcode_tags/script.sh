#!/bin/bash

set -eo pipefail

#  Sam tags added by BD Rhapsody Pipeline
#  From: https://www.bd.com/documents/guides/user-guides/GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf
# 
# =========================================================================================
# |    | Definition                                                                       |
# =========================================================================================
# | CB | A number between 1 and 96 3 (884,736) representing a unique cell label sequence  |
# |    | (CB = 0 when no cell label sequence is detected)                                 |
# -----------------------------------------------------------------------------------------
# | MR | Raw molecular identifier sequence                                                |
# -----------------------------------------------------------------------------------------
# | MA | RSEC-adjusted molecular identifier sequence. If not a true cell, the raw UMI is  |
# |    | repeated in this tag.                                                            |
# -----------------------------------------------------------------------------------------
# | PT | T if a poly(T) tail was found in the expected position on R1, or F if poly(T)    |
# |    | was not found                                                                    |
# -----------------------------------------------------------------------------------------
# | CN | Indicates if a sequence is derived from a putative cell, as determined by the    |
# |    | cell label filtering algorithm (T: putative cell; x: invalid cell label or noise |
# |    | cell) Note: You can distinguish between an invalid cell label and a noise cell   |
# |    | with the CB tag (invalid cell labels are 0).                                     |
# -----------------------------------------------------------------------------------------
# | ST | The value is 1-12, indicating the Sample Tag of the called putative cell, or M   |
# |    | for multiplet, or x for undetermined.                                            |
# =========================================================================================


# SAM tags added by 10X
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
# =========================================================================================
# |    | Definition                                                                       |
# =========================================================================================
# | CB | Chromium cellular barcode sequence that is error-corrected and confirmed against |
# |    | a list of known-good barcode sequences. For multiplex Fixed RNA Profiling, the   |
# |    | cellular barcode is a combination of the 10x GEM Barcode and Probe Barcode       |
# |    | sequences.                                                                       |
# -----------------------------------------------------------------------------------------
# | CR | Chromium cellular barcode sequence as reported by the sequencer. For multiplex   |
# |    | Fixed RNA Profiling, the cellular barcode is a combination of the 10x GEM        |
# |    | Barcode and Probe Barcode sequences.                                             |
# -----------------------------------------------------------------------------------------
# | CY | Chromium cellular barcode read quality. For multiplex Fixed RNA Profiling, the   |
# |    | cellular barcode is a combination of the 10x GEM Barcode and Probe Barcode       |
# |    | sequences. Phred scores as reported by sequencer.                                |
# -----------------------------------------------------------------------------------------
# | UB | Chromium molecular barcode sequence that is error-corrected among other          |
# |    | molecular barcodes with the same cellular barcode and gene alignment.            |
# -----------------------------------------------------------------------------------------
# | UR | Chromium molecular barcode sequence as reported by the sequencer.                |
# -----------------------------------------------------------------------------------------
# | UY | Chromium molecular barcode read quality. Phred scores as reported by sequencer.  |
# -----------------------------------------------------------------------------------------
# | TR | Trimmed sequence. For the Single Cell 3' v1 chemistry, this is trailing sequence |
# |    | following the UMI on Read 2. For the Single Cell 3' v2 chemistry, this is        |
# |    | trailing sequence following the cell and molecular barcodes on Read 1.           |
# =========================================================================================

extra_params=()

if [ "$par_bam" == "true" ]; then
  extra_params+=("--bam")
fi

cat \
    <(samtools view -SH "$par_input") \
    <(samtools view "$par_input" |  grep "MA:Z:*"  | sed  "s/MA:Z:/UB:Z:/" ) | \
samtools view -Sh "${extra_params[@]}" -@"$par_threads" - > "$par_output"