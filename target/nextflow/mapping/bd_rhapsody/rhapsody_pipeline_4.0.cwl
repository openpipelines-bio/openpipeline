#!/usr/bin/env cwl-runner
{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 64000
                }
            ],
            "baseCommand": "ATAC_Cell_by_Peak.R",
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--fragments"
                    },
                    "secondaryFiles": {
                        "pattern": ".tbi",
                        "required": null
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/Fragments"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/GTF"
                },
                {
                    "type": "boolean",
                    "inputBinding": {
                        "prefix": "--map-gene-name-id"
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/Map_Gene_Name_ID"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--peaks"
                    },
                    "secondaryFiles": {
                        "pattern": ".tbi",
                        "required": null
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/Peaks"
                },
                {
                    "type": "File",
                    "loadContents": true,
                    "id": "#ATAC_Cell_by_Peak.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--transposase-sites"
                    },
                    "secondaryFiles": {
                        "pattern": ".tbi",
                        "required": null
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/Transposase_Sites"
                }
            ],
            "arguments": [
                {
                    "prefix": "--base-name",
                    "valueFrom": "$(JSON.parse(inputs.Run_Metadata.contents).Run_Base_Name)"
                }
            ],
            "id": "#ATAC_Cell_by_Peak.cwl",
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_ATAC_Cell_Calling_Data.csv"
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/ATAC_Cell_Calling_Data"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*atac.gene_info.json"
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/Gene_Name_ID_Json"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Initial_Seurat.rds"
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/Initial_Seurat_RDS"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Total_Fragment_Metrics.json"
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/Total_Fragment_Metrics"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#ATAC_Cell_by_Peak.cwl/output"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 48000
                }
            ],
            "baseCommand": "ATAC_Compile_Results.R",
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--biop-putative-data-table"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Biop_putative_data_table"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--cell-order"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Cell_Order"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--cell-order-subsampled"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Cell_Order_Subsampled"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--fragments"
                    },
                    "secondaryFiles": {
                        "pattern": ".tbi",
                        "required": null
                    },
                    "id": "#ATAC_Compile_Results.cwl/Fragments"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--feature-name-json"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Gene_Name_ID_Json"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--initial-seurat"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Initial_Seurat_RDS"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--input-metrics-tar"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Input_Metrics_tar"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--genome-contig-lengths"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Reference_Genome_Contig_Lengths"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "prefix": "--genome-size"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Reference_Genome_Size"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--total-fragment-metrics"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Total_Fragment_Metrics"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--unified-metrics"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Unified_Metrics"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "*_MEX.zip"
                    },
                    "id": "#ATAC_Compile_Results.cwl/ATAC_Data_Tables"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "*/*_coordinates.csv"
                    },
                    "id": "#ATAC_Compile_Results.cwl/ATAC_Dim_Reduction_Coord"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_ATAC_Nucleosome_Signal.csv"
                    },
                    "id": "#ATAC_Compile_Results.cwl/ATAC_Nucleosome_Signal"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_ATAC_Seurat.rds"
                    },
                    "id": "#ATAC_Compile_Results.cwl/ATAC_Seurat_RDS"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_ATAC_TSS_Enrichment.csv"
                    },
                    "id": "#ATAC_Compile_Results.cwl/ATAC_TSS_Enrichment"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_ATAC_Transposase_Sites_In_Genome.csv"
                    },
                    "id": "#ATAC_Compile_Results.cwl/ATAC_Transposase_Site_Counts"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*cell_type_experimental.csv"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Cell_Type_Predictions"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "metrics-files.tar.gz"
                    },
                    "id": "#ATAC_Compile_Results.cwl/Metrics_tar"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "mist_atac_compile_results.log"
                    },
                    "id": "#ATAC_Compile_Results.cwl/output"
                }
            ],
            "id": "#ATAC_Compile_Results.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 48000
                }
            ],
            "baseCommand": [
                "mist_atac_motif_analysis.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--atac-cell-by-peak"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/Atac_Datatables"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/Maximum_Threads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--reference-fasta"
                    },
                    "secondaryFiles": [
                        {
                            "pattern": ".fai",
                            "required": null
                        },
                        {
                            "pattern": ".gzi",
                            "required": null
                        }
                    ],
                    "id": "#ATAC_Motif_Analysis.cwl/Reference_FASTA"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--pfm-file"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/Transcription_Factor_Motif_PFM"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "ATAC-Background-Peaks.csv.gz"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/ATAC_Background_Peaks"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "*_MEX.zip"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/ATAC_Motif_Data_Tables"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_ATAC_Motif_Positions.tsv.gz"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/ATAC_Motif_Positions"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "ATAC-Peak-Seq.csv.gz"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/ATAC_Peak_Seq"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#ATAC_Motif_Analysis.cwl/output"
                }
            ],
            "id": "#ATAC_Motif_Analysis.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "coresMin": "${ if (inputs.Assay == 'ATAC') { return 2; } else { return 4; } }",
                    "ramMin": "${ if (inputs.Assay == 'ATAC') { return 4000; } else { return 32000; } }"
                }
            ],
            "baseCommand": [
                "mist_add_to_bam.py"
            ],
            "inputs": [
                {
                    "type": "string",
                    "inputBinding": {
                        "prefix": "--assay"
                    },
                    "id": "#AddtoBam.cwl/Assay"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--input-bam"
                    },
                    "id": "#AddtoBam.cwl/Bam"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--cell-order-json"
                    },
                    "id": "#AddtoBam.cwl/Cell_Order"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "prefix": "--corrected-mols-list",
                        "itemSeparator": ","
                    },
                    "id": "#AddtoBam.cwl/Corrected_Mols"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#AddtoBam.cwl/Generate_Bam"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata-json"
                    },
                    "id": "#AddtoBam.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--sample-tag-calls"
                    },
                    "id": "#AddtoBam.cwl/SampleTag_Calls"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--target-gene-mapping"
                    },
                    "id": "#AddtoBam.cwl/Target_Gene_Mapping"
                }
            ],
            "arguments": [
                {
                    "prefix": "--bamIO-threads",
                    "valueFrom": "$(runtime.cores)"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.bam"
                    },
                    "id": "#AddtoBam.cwl/Annotated_Bam"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#AddtoBam.cwl/output"
                }
            ],
            "id": "#AddtoBam.cwl"
        },
        {
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "doc": "AlignmentAnalysis stage of the Rain pipeline annotates aligned reads and collects a myriad of metrics on the aligned reads. Additional annotation is performed to the reads\n",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "coresMin": 8,
                    "ramMin": 24000
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "AlignmentAnalysisAndCountCB.sh"
            ],
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--assay"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Assay"
                },
                {
                    "inputBinding": {
                        "prefix": "--exclude-intronic-reads"
                    },
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Exclude_Intronic_Reads"
                },
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Extra_Seqs"
                },
                {
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/GTF"
                },
                {
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Maximum_Threads"
                },
                {
                    "inputBinding": {
                        "prefix": "--multiple-alignment-gene-selection"
                    },
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Multiple_Alignment_Gene_Selection"
                },
                {
                    "inputBinding": {
                        "prefix": "--multiple-alignment-handling"
                    },
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Multiple_Alignment_Handling"
                },
                {
                    "inputBinding": {
                        "prefix": "--r2-bam",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#AlignmentAnalysis.cwl/R2_BAM"
                },
                {
                    "inputBinding": {
                        "prefix": "--quality-metrics"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/ReadQualityMetrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#AlignmentAnalysis.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--transcript-length"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Transcript_Length"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*.annotated.*.bam"
                    },
                    "id": "#AlignmentAnalysis.cwl/Annotated_Bam_Files"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*gtf.gene_info.json"
                    },
                    "id": "#AlignmentAnalysis.cwl/Gene_Name_ID_Json"
                },
                {
                    "outputBinding": {
                        "glob": "*logs.tar.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignmentAnalysis.cwl/Logs"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_SeqMetrics.csv"
                    },
                    "id": "#AlignmentAnalysis.cwl/Seq_Metrics"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*Sorted_Valid_Reads.csv.*"
                    },
                    "id": "#AlignmentAnalysis.cwl/Sorted_Valid_Reads_CSV"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "count_estimates.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 30000; } return parseInt(JSON.parse(self[0].contents).num_bioproducts); }"
                    },
                    "id": "#AlignmentAnalysis.cwl/num_bioproducts"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "count_estimates.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 10000; } return parseInt(JSON.parse(self[0].contents).num_cell_estimate); }"
                    },
                    "id": "#AlignmentAnalysis.cwl/num_cell_estimate"
                },
                {
                    "type": "long",
                    "outputBinding": {
                        "glob": "num_vdj_reads.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 0; } return parseInt(JSON.parse(self[0].contents).BCR); }"
                    },
                    "id": "#AlignmentAnalysis.cwl/num_valid_ig_reads"
                },
                {
                    "type": "long",
                    "outputBinding": {
                        "glob": "num_vdj_reads.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 0; } return parseInt(JSON.parse(self[0].contents).TCR); }"
                    },
                    "id": "#AlignmentAnalysis.cwl/num_valid_tcr_reads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_BCR_Valid_Reads.fastq.gz"
                    },
                    "id": "#AlignmentAnalysis.cwl/validIgReads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_TCR_Valid_Reads.fastq.gz"
                    },
                    "id": "#AlignmentAnalysis.cwl/validTcrReads"
                }
            ],
            "id": "#AlignmentAnalysis.cwl"
        },
        {
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_annotate_molecules.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--umi-option"
                    },
                    "id": "#AnnotateMolecules.cwl/AbSeq_UMI"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#AnnotateMolecules.cwl/Assay"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#AnnotateMolecules.cwl/Run_Metadata"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--valid-annot"
                    },
                    "id": "#AnnotateMolecules.cwl/Valids"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_bioproduct_stats*.json"
                    },
                    "id": "#AnnotateMolecules.cwl/Bioproduct_Stats_List"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_CellBiopSummary.csv.*"
                    },
                    "id": "#AnnotateMolecules.cwl/Cell_Biop_Summary_List"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_Annotation_Molecule_corrected.csv.*"
                    },
                    "id": "#AnnotateMolecules.cwl/Corrected_Mols_List"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "stats.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).num_molecules_corrected)\n"
                    },
                    "id": "#AnnotateMolecules.cwl/Num_Molecules_Corrected"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "stats.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).total_molecules)\n"
                    },
                    "id": "#AnnotateMolecules.cwl/Total_Molecules"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#AnnotateMolecules.cwl/output"
                }
            ],
            "id": "#AnnotateMolecules.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "id": "#Assay_Settings.cwl/AbSeq_Reference_Present"
                },
                {
                    "type": "int",
                    "id": "#Assay_Settings.cwl/Cell_Calling_Data"
                },
                {
                    "type": "boolean",
                    "id": "#Assay_Settings.cwl/Reads_ATAC_Present"
                },
                {
                    "type": "boolean",
                    "id": "#Assay_Settings.cwl/Reads_RNA_Present"
                },
                {
                    "type": "boolean",
                    "id": "#Assay_Settings.cwl/Reference_Archive_Present"
                },
                {
                    "type": "boolean",
                    "id": "#Assay_Settings.cwl/Targeted_Reference_Present"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Assay_Settings.cwl/Assay_ATAC"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Assay_Settings.cwl/Assay_RNA"
                }
            ],
            "expression": "${\n  var assay_rna = null;\n  var assay_atac = null;\n\n  if (!inputs.Reads_ATAC_Present && !inputs.Reads_RNA_Present)\n  {\n    throw new Error('Invalid pipeline inputs: Please provide Reads for at least 1 of RNA or ATAC analysis.')\n  }\n  if (inputs.Targeted_Reference_Present && inputs.Reference_Archive_Present) {\n    throw new Error('Invalid pipeline inputs: Do not provide both Targeted Reference and Reference Archive.')\n  }\n  if (!inputs.Targeted_Reference_Present && !inputs.AbSeq_Reference_Present && !inputs.Reference_Archive_Present && inputs.Cell_Calling_Data != 6) {\n    throw new Error('Invalid pipeline inputs: Please provide either a Reference Archive or a Targeted Reference or an AbSeq Reference. If you are analyzing VDJ data only, then please choose VDJ only cell calling algorithm.')\n  }\n\n  if ( inputs.Reads_ATAC_Present )\n  {\n    assay_atac = 'ATAC'\n  }\n\n  if (inputs.Reads_RNA_Present && inputs.Reference_Archive_Present && inputs.Cell_Calling_Data != 6) {\n    assay_rna = 'WTA'\n  }\n  else if (inputs.Reads_RNA_Present && (inputs.Targeted_Reference_Present || inputs.AbSeq_Reference_Present || inputs.Cell_Calling_Data == 6)) {\n    assay_rna = 'Targeted'\n  }\n\n\n  return ({Assay_RNA: assay_rna, Assay_ATAC: assay_atac})\n}",
            "id": "#Assay_Settings.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#BamSettings.cwl/_Generate_Bam"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#BamSettings.cwl/Generate_Bam"
                }
            ],
            "expression": "${\n  // the create bam flag defaults to false\n  var generateBam = false;\n  // the user can set this flag to true, to enable creation of the bam file.\n  if (inputs._Generate_Bam != null) {\n    generateBam = inputs._Generate_Bam;\n  }\n  return ({\n    Generate_Bam: generateBam,\n  });\n}",
            "id": "#BamSettings.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#BundleLogs.cwl/log_files"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "id": "#BundleLogs.cwl/logs_dir"
                }
            ],
            "expression": "${\n  /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */\n  function uuid() {\n    var uuid = \"\", i, random;\n    for (i = 0; i < 32; i++) {\n      random = Math.random() * 16 | 0;\n      if (i == 8 || i == 12 || i == 16 || i == 20) {\n        uuid += \"-\";\n      }\n      uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);\n    }\n    return uuid;\n  }\n  var listing = [];\n  for (var i = 0; i < inputs.log_files.length; i++) {\n    var log_file = inputs.log_files[i];\n    /*\n      Checking here for null in case a Node was skipped because of conditional execution.\n      For e.g. Generate_Bam is used to skip the AddToBam, MergeBam and IndexBam nodes\n    */\n    if (log_file != null) {\n      log_file.basename = uuid() + \"-\" + log_file.basename;\n      listing.push(log_file);\n    }\n  }\n  return ({\n    logs_dir: {\n      class: \"Directory\",\n      basename: \"Logs\",\n      listing: listing\n    }\n  });\n}",
            "id": "#BundleLogs.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "This step computes gene body coverage, using only exonic reads from the bam file and exons from the GTF file.  It will output the result in a json file that has 100 values. These 100 value represents the mean normalized coverage  in 0-100 percentile of a gene body length.\n",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "coresMin": 8
                }
            ],
            "baseCommand": [
                "mist_calc_gene_body_coverage.py"
            ],
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--assay"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#CalcGeneBodyCoverage.cwl/Assay"
                },
                {
                    "inputBinding": {
                        "prefix": "--bam-file-list",
                        "itemSeparator": ",",
                        "separate": true
                    },
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#CalcGeneBodyCoverage.cwl/Bam"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "id": "#CalcGeneBodyCoverage.cwl/Gtf"
                },
                {
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CalcGeneBodyCoverage.cwl/Maximum_Threads"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "prefix": "--run-name"
                    },
                    "id": "#CalcGeneBodyCoverage.cwl/Run_Base_Name"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_gene_body_coverage.json"
                    },
                    "id": "#CalcGeneBodyCoverage.cwl/Gene_Body_Coverage_Json"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#CalcGeneBodyCoverage.cwl/output"
                }
            ],
            "id": "#CalcGeneBodyCoverage.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "doc": "This step computes sequencing saturation at various downsampled mean useful reads per cell values and outputs the necessary information to construct a sequencing saturation curve\n",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": "${\n  const BASE_MEM_FOR_NODE = 2 * 1024 * 1024 * 1024; // Base : 2G memory\n  const BYTES_PER_ENTRY = 128;\n  const totalLines = inputs.Num_Mols_Corrected || 0; // Then scale based on the # rsec corrected mols - used to contruct a look-up table\n  const requiredBytes = totalLines * BYTES_PER_ENTRY + BASE_MEM_FOR_NODE;\n  return Math.ceil(requiredBytes / (1024 * 1024));\n}\n",
                    "coresMin": "${\n  const n = 8;\n  if (inputs.Threads) {\n    return Math.min(inputs.Threads, n);\n  } else {\n    return n;\n  }\n}\n"
                }
            ],
            "baseCommand": [
                "mist_calc_seq_saturation.py"
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "id": "#CalcSeqSaturation.cwl/AbSeq_Reference_Present"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "prefix": "--corrected-mols-list",
                        "itemSeparator": ","
                    },
                    "id": "#CalcSeqSaturation.cwl/Corrected_Mols"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--metrics-tar"
                    },
                    "id": "#CalcSeqSaturation.cwl/Metrics_Tar"
                },
                {
                    "type": "int",
                    "id": "#CalcSeqSaturation.cwl/Num_Mols_Corrected"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--putative-cells-json-fp"
                    },
                    "id": "#CalcSeqSaturation.cwl/Putative_Cells"
                },
                {
                    "type": "boolean",
                    "id": "#CalcSeqSaturation.cwl/Reference_Archive_Present"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "prefix": "--run-name"
                    },
                    "id": "#CalcSeqSaturation.cwl/Run_Base_Name"
                },
                {
                    "type": "boolean",
                    "id": "#CalcSeqSaturation.cwl/Targeted_Reference_Present"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#CalcSeqSaturation.cwl/Threads"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "prefix": "--valid-reads-csv-list",
                        "itemSeparator": ","
                    },
                    "id": "#CalcSeqSaturation.cwl/Valid_Reads_CSV"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_seq_saturation_plot.json"
                    },
                    "id": "#CalcSeqSaturation.cwl/Seq_Saturation_Plot_Json"
                }
            ],
            "id": "#CalcSeqSaturation.cwl"
        },
        {
            "requirements": [
                {
                    "listing": [
                        "$(inputs.Reference_Archive)"
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_check_references.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--abseq-reference"
                    },
                    "id": "#CheckReference.cwl/AbSeq_Reference"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--assay"
                    },
                    "id": "#CheckReference.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "id": "#CheckReference.cwl/Cell_Calling_Data"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--predefined-peaks"
                    },
                    "id": "#CheckReference.cwl/Predefined_ATAC_Peaks"
                },
                {
                    "type": [
                        "null",
                        "File",
                        "Directory"
                    ],
                    "inputBinding": {
                        "prefix": "--reference-archive"
                    },
                    "id": "#CheckReference.cwl/Reference_Archive"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--sample-tags-version"
                    },
                    "id": "#CheckReference.cwl/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--supplemental-reference"
                    },
                    "id": "#CheckReference.cwl/Supplemental_Reference"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--targeted-reference"
                    },
                    "id": "#CheckReference.cwl/Targeted_Reference"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-version"
                    },
                    "id": "#CheckReference.cwl/VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_ATAC_Peaks.bed.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "secondaryFiles": [
                        {
                            "pattern": ".tbi",
                            "required": null
                        }
                    ],
                    "id": "#CheckReference.cwl/Checked_Predefined_Peaks"
                },
                {
                    "outputBinding": {
                        "glob": "*_ATAC_Peaks.bed.gz.tbi"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Checked_Predefined_Peaks_Index"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "combined_extra_seq.fasta"
                    },
                    "id": "#CheckReference.cwl/Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "full-gene-list.json"
                    },
                    "id": "#CheckReference.cwl/Full_Genes"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "BD_Rhapsody_Reference_Files/*.gtf"
                    },
                    "id": "#CheckReference.cwl/GTF"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "BD_Rhapsody_Reference_Files"
                    },
                    "id": "#CheckReference.cwl/Index"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "BD_Rhapsody_Reference_Files/*.fasta.gz"
                    },
                    "secondaryFiles": [
                        {
                            "pattern": ".fai",
                            "required": null
                        },
                        {
                            "pattern": ".gzi",
                            "required": null
                        }
                    ],
                    "id": "#CheckReference.cwl/Reference_FASTA"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "target-gene.json"
                    },
                    "id": "#CheckReference.cwl/Target_Gene_Mapping"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "transcript_length.json"
                    },
                    "id": "#CheckReference.cwl/Transcript_Length"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "BD_Rhapsody_Reference_Files/*.pfm"
                    },
                    "id": "#CheckReference.cwl/Transcription_Factor_Motif_PFM"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#CheckReference.cwl/output"
                }
            ],
            "id": "#CheckReference.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 32000
                }
            ],
            "baseCommand": [
                "mist_generate_H5MU.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-background-peaks"
                    },
                    "id": "#GenerateH5MU.cwl/Atac_Background_Peaks"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--atac-cell-by-peak"
                    },
                    "id": "#GenerateH5MU.cwl/Atac_Datatables"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--atac-metrics"
                    },
                    "id": "#GenerateH5MU.cwl/Atac_Metrics"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-motif-positions"
                    },
                    "id": "#GenerateH5MU.cwl/Atac_Motif_Positions"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-nucleosome-signal"
                    },
                    "id": "#GenerateH5MU.cwl/Atac_Nucleosome_Signal"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-peak-seq"
                    },
                    "id": "#GenerateH5MU.cwl/Atac_Peak_Seq"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-tss-enrichment"
                    },
                    "id": "#GenerateH5MU.cwl/Atac_TSS_Enrichment"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-transposase-site-counts"
                    },
                    "id": "#GenerateH5MU.cwl/Atac_Transposase_Site_Counts"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--bioproduct-stats"
                    },
                    "id": "#GenerateH5MU.cwl/Bioproduct_Stats"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "prefix": "--cell-type-experimental",
                        "itemSeparator": ","
                    },
                    "id": "#GenerateH5MU.cwl/Cell_Type_Experimental"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--data-tables"
                    },
                    "id": "#GenerateH5MU.cwl/Data_Tables"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--coordinates-file-list"
                    },
                    "id": "#GenerateH5MU.cwl/Dim_Reduction_Coord"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "id": "#GenerateH5MU.cwl/GTF"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--metrics-summary"
                    },
                    "id": "#GenerateH5MU.cwl/Metrics_Summary"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--peak_annotation"
                    },
                    "id": "#GenerateH5MU.cwl/Peak_Annotation"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--pipeline-report-html"
                    },
                    "id": "#GenerateH5MU.cwl/Pipeline_Report_Html"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--protein-aggregates-experimental"
                    },
                    "id": "#GenerateH5MU.cwl/Protein_Aggregates_Experimental"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--putative-cells-origin"
                    },
                    "id": "#GenerateH5MU.cwl/Putative_Cells_Origin"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--genome-contig-lengths"
                    },
                    "id": "#GenerateH5MU.cwl/Reference_Genome_Contig_Lengths"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#GenerateH5MU.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--sample-tag-files"
                    },
                    "id": "#GenerateH5MU.cwl/SampleTag_CSVs"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--sample-tag-calls"
                    },
                    "id": "#GenerateH5MU.cwl/SampleTag_Calls"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--pfm-file"
                    },
                    "id": "#GenerateH5MU.cwl/Transcription_Factor_Motif_PFM"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-anndata"
                    },
                    "id": "#GenerateH5MU.cwl/VDJ_AnnData"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-per-cell"
                    },
                    "id": "#GenerateH5MU.cwl/VDJ_Per_Cell"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.cellismo"
                    },
                    "id": "#GenerateH5MU.cwl/H5MU"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Pipeline_Report.html"
                    },
                    "id": "#GenerateH5MU.cwl/Pipeline_Report_HTML_Inj"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_biotypeCounts.csv"
                    },
                    "id": "#GenerateH5MU.cwl/biotypeCounts"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#GenerateH5MU.cwl/output"
                }
            ],
            "id": "#GenerateH5MU.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 32000
                }
            ],
            "baseCommand": "GenerateSeurat.R",
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-seurat-rds"
                    },
                    "id": "#GenerateSeurat.cwl/ATAC_Seurat"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-background-peaks"
                    },
                    "id": "#GenerateSeurat.cwl/Atac_Background_Peaks"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--atac-data-tables"
                    },
                    "id": "#GenerateSeurat.cwl/Atac_Datatables"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--atac-metrics"
                    },
                    "id": "#GenerateSeurat.cwl/Atac_Metrics"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--atac-motif-data-table"
                    },
                    "id": "#GenerateSeurat.cwl/Atac_Motif_Datatables"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-motif-positions"
                    },
                    "id": "#GenerateSeurat.cwl/Atac_Motif_Positions"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-peak-seq"
                    },
                    "id": "#GenerateSeurat.cwl/Atac_Peak_Seq"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--bioproduct-stats"
                    },
                    "id": "#GenerateSeurat.cwl/Bioproduct_Stats"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "prefix": "--cell-type-experimental",
                        "itemSeparator": ","
                    },
                    "id": "#GenerateSeurat.cwl/Cell_Type_Experimental"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--data-tables"
                    },
                    "id": "#GenerateSeurat.cwl/Data_Tables"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--coordinates-file-list"
                    },
                    "id": "#GenerateSeurat.cwl/Dim_Reduction_Coord"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--gtf-file"
                    },
                    "id": "#GenerateSeurat.cwl/GTF"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--metrics-summary"
                    },
                    "id": "#GenerateSeurat.cwl/Metrics_Summary"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--peak-annotation"
                    },
                    "id": "#GenerateSeurat.cwl/Peak_Annotation"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--protein-aggregates-experimental"
                    },
                    "id": "#GenerateSeurat.cwl/Protein_Aggregates_Experimental"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--putative-cells-origin"
                    },
                    "id": "#GenerateSeurat.cwl/Putative_Cells_Origin"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--genome-contig-lengths"
                    },
                    "id": "#GenerateSeurat.cwl/Reference_Genome_Contig_Lengths"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--metadata-file"
                    },
                    "id": "#GenerateSeurat.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--sample-tag-csvs"
                    },
                    "id": "#GenerateSeurat.cwl/SampleTag_CSVs"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--sample-tag-calls"
                    },
                    "id": "#GenerateSeurat.cwl/SampleTag_Calls"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--pfm-file"
                    },
                    "id": "#GenerateSeurat.cwl/Transcription_Factor_Motif_PFM"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-per-cell"
                    },
                    "id": "#GenerateSeurat.cwl/VDJ_Per_Cell"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.rds"
                    },
                    "id": "#GenerateSeurat.cwl/SeuratRDS"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#GenerateSeurat.cwl/output"
                }
            ],
            "id": "#GenerateSeurat.cwl"
        },
        {
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 48000
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_get_datatables.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--atac-cell-calling-data-file"
                    },
                    "id": "#GetDataTable.cwl/ATAC_Cell_Calling_Input"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--bioproduct-stats-list"
                    },
                    "id": "#GetDataTable.cwl/Bioproduct_Stats_List"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--cell-biop-summary-list"
                    },
                    "id": "#GetDataTable.cwl/Cell_Biop_Summary_List"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--full-gene-list"
                    },
                    "id": "#GetDataTable.cwl/Full_Genes"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--feature-name-json"
                    },
                    "id": "#GetDataTable.cwl/Gene_Name_ID_Json"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#GetDataTable.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--seq-metrics"
                    },
                    "id": "#GetDataTable.cwl/Seq_Metrics"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "id": "#GetDataTable.cwl/Total_Molecules"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#GetDataTable.cwl/num_bioproducts"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#GetDataTable.cwl/num_cell_estimate"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_RSEC_MolsPerCell_MEX.zip"
                    },
                    "id": "#GetDataTable.cwl/Biop_putative_data_table"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Bioproduct_Stats.csv"
                    },
                    "id": "#GetDataTable.cwl/Bioproduct_Stats"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "cell_order.json"
                    },
                    "id": "#GetDataTable.cwl/Cell_Order"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "cell_order_subsampled.json"
                    },
                    "id": "#GetDataTable.cwl/Cell_Order_Subsampled"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*cell_type_experimental.csv"
                    },
                    "id": "#GetDataTable.cwl/Cell_Type_Predictions"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "*_MEX.zip"
                    },
                    "id": "#GetDataTable.cwl/Data_Tables"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "*_coordinates.csv"
                    },
                    "id": "#GetDataTable.cwl/Dim_Reduction_Coord"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "metrics-files.tar.gz"
                    },
                    "id": "#GetDataTable.cwl/Metrics_tar"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Protein_Agg/*_Protein_Aggregates_Experimental.csv"
                    },
                    "id": "#GetDataTable.cwl/Protein_Aggregates_Experimental"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*_Putative_Cells_Origin.csv"
                    },
                    "id": "#GetDataTable.cwl/Putative_Cells_Origin"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "SampleTag/*csv"
                    },
                    "id": "#GetDataTable.cwl/SampleTag_CSVs"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "SampleTag/*_Sample_Tag_Calls.csv"
                    },
                    "id": "#GetDataTable.cwl/SampleTag_Calls"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "SampleTagArchives/*zip"
                    },
                    "id": "#GetDataTable.cwl/SampleTag_perTagZips"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#GetDataTable.cwl/output"
                }
            ],
            "id": "#GetDataTable.cwl"
        },
        {
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "doc": "GetMachineResources gets available resources (current only the number of cpus) on the machine that is running the local deployment of the MIST pipeline.\n",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_machine_resources.py"
            ],
            "inputs": [],
            "outputs": [
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "resources.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).total_cpus_avail)"
                    },
                    "id": "#GetMachineResources.cwl/Total_CPUs_Avail"
                }
            ],
            "id": "#GetMachineResources.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/_AbSeq_UMI"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/_MinChunkSize"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#InternalSettings.cwl/_NumRecordsPerSplit"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/_Subsample_Sample_Tags"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/_Target_analysis"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/_VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/_VDJ_VGene_Evalue"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/AbSeq_UMI"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/MinChunkSize"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#InternalSettings.cwl/NumRecordsPerSplit"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/Subsample_Sample_Tags"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Target_analysis"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/VDJ_VGene_Evalue"
                }
            ],
            "expression": "${\n  var internalInputs = [\n    '_AbSeq_UMI',\n    '_MinChunkSize',\n    '_NumRecordsPerSplit',\n    '_Target_analysis',\n    '_Subsample_Sample_Tags',\n    '_VDJ_VGene_Evalue',\n    '_VDJ_JGene_Evalue',\n  ];\n  var internalOutputs = {}\n  for (var i = 0; i < internalInputs.length; i++) {\n    var internalInput = internalInputs[i];\n    var internalOutput = internalInput.slice(1); // remove leading underscore\n    if (inputs.hasOwnProperty(internalInput)) {\n      internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output\n    } else {\n      internalOutputs[internalOutput] = null; // if input not specified, provide a null\n    }\n  }\n  return internalOutputs;\n}",
            "id": "#InternalSettings.cwl"
        },
        {
            "class": "Workflow",
            "label": "BD Rhapsody\u2122 Sequence Analysis Pipeline",
            "doc": "The BD Rhapsody\u2122 assays are used to create sequencing libraries from single cell transcriptomes.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files and a reference file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "AbSeq Reference",
                    "id": "#main/AbSeq_Reference"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/AbSeq_UMI"
                },
                {
                    "label": "Cell Calling ATAC Algorithm",
                    "doc": "Specify the ATAC algorithm to be used for ATAC putative cell calling. The Basic algorithm is the default.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Cell_Calling_ATAC_Algorithm/Cell_Calling_ATAC_Algorithm",
                            "symbols": [
                                "#main/Cell_Calling_ATAC_Algorithm/Cell_Calling_ATAC_Algorithm/Basic",
                                "#main/Cell_Calling_ATAC_Algorithm/Cell_Calling_ATAC_Algorithm/Refined"
                            ]
                        }
                    ],
                    "id": "#main/Cell_Calling_ATAC_Algorithm"
                },
                {
                    "label": "Cell Calling Bioproduct Algorithm",
                    "doc": "Specify the bioproduct algorithm to be used for mRNA/AbSeq putative cell calling. The Basic algorithm is the default.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Cell_Calling_Bioproduct_Algorithm/Cell_Calling_Bioproduct_Algorithm",
                            "symbols": [
                                "#main/Cell_Calling_Bioproduct_Algorithm/Cell_Calling_Bioproduct_Algorithm/Basic",
                                "#main/Cell_Calling_Bioproduct_Algorithm/Cell_Calling_Bioproduct_Algorithm/Refined"
                            ]
                        }
                    ],
                    "id": "#main/Cell_Calling_Bioproduct_Algorithm"
                },
                {
                    "label": "Cell Calling Data",
                    "doc": "Specify the data to be used for putative cell calling.\nThe default data for putative cell calling will be determined the following way:\n  - If mRNA and ATAC Reads exist, mRNA_and_ATAC is the default.\n  - If only ATAC Reads exist, ATAC is the default.\n  - Otherwise, mRNA is the default.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Cell_Calling_Data/Cell_Calling_Data",
                            "symbols": [
                                "#main/Cell_Calling_Data/Cell_Calling_Data/mRNA",
                                "#main/Cell_Calling_Data/Cell_Calling_Data/AbSeq",
                                "#main/Cell_Calling_Data/Cell_Calling_Data/ATAC",
                                "#main/Cell_Calling_Data/Cell_Calling_Data/mRNA_and_ATAC",
                                "#main/Cell_Calling_Data/Cell_Calling_Data/VDJ"
                            ]
                        }
                    ],
                    "id": "#main/Cell_Calling_Data"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Custom STAR Params",
                    "doc": "Allows you to specify custom STAR aligner mapping parameters. Only the mapping parameters you provide here will be used with STAR, meaning that you must provide the complete list of parameters that you want to take effect. For reference, the parameters used by default in the pipeline are:\n\n    1. Short Reads: --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMultimapScoreRange 0 --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAA --seedSearchStartLmax 50 --outFilterMatchNmin 25 --limitOutSJcollapsed 2000000\n    2. Long Reads: Same options as short reads + --seedPerReadNmax 10000\n\n\nExample input: --alignIntronMax 500000 --outFilterScoreMinOverLread 0 --limitOutSJcollapsed 2000000\n\nImportant:\n  1. This applies to fastqs provided in the Reads user input\n  2. Please do not specify any non-mapping related params like: --runThreadN, --genomeDir --outSAMtype, etc.\n  3. Please only use params supported by STAR version 2.7.10b\n",
                    "id": "#main/Custom_STAR_Params"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Custom bwa-mem2 Params",
                    "doc": "Allows you to specify custom bwa-mem2 mapping parameters. Only the mapping parameters you provide here will be used with bwa-mem2, meaning that you must provide the complete list of parameters that you want to take effect. The pipeline uses program default mapping parameters.\n\nExample input: -k 15 -w 200 -r 2\n\nImportant:\n    1. This applies to fastqs provided in the Reads_ATAC user input\n    2. Please do not specify any non-mapping related params like: -C, -t, etc.\n    3. Please only use params supported by bwa-mem2 version 2.2.1\n",
                    "id": "#main/Custom_bwa_mem2_Params"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Exact Cell Count",
                    "doc": "Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count",
                    "id": "#main/Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Exclude Intronic Reads",
                    "doc": "By default, reads aligned to exons and introns are considered and represented in molecule counts. Including intronic reads may increase sensitivity, resulting in an increase in molecule counts and the number of genes per cell for both cellular and nuclei samples. Intronic reads may indicate unspliced mRNAs and are also useful, for example, in the study of nuclei and RNA velocity. When set to true, intronic reads will be excluded.",
                    "id": "#main/Exclude_Intronic_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Expected Cell Count",
                    "doc": "Optional.  Guide the basic putative cell calling algorithm by providing an estimate of the number of cells expected.  Usually this can be the number of cells loaded into the Rhapsody cartridge.  If there are multiple inflection points on the second derivative cumulative curve, this will ensure the one selected is near the expected.",
                    "id": "#main/Expected_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Generate Bam Output",
                    "doc": "Default: false.  A Bam read alignment file contains reads from all the input libraries, but creating it can consume a lot of compute and disk resources. By setting this field to true, the Bam file will be created. This option is shared for both Bioproduct and ATAC libraries.\n",
                    "id": "#main/Generate_Bam"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Long Reads (>=650bp)",
                    "doc": "By default, we detect if there are any reads longer than 650bp and then flag QualCLAlign to use STARlong instead of STAR. This flag can be explicitly set if it is known in advance that there are reads longer than 650bp.\n",
                    "id": "#main/Long_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Maximum Number of Threads",
                    "doc": "The maximum number of threads to use in the pipeline. By default, all available cores are used.",
                    "id": "#main/Maximum_Threads"
                },
                {
                    "label": "Multiple Alignment Gene Selection",
                    "doc": "Specify the feature selection strategy for multi-aligned reads.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Multiple_Alignment_Gene_Selection/Multiple_Alignment_Gene_Selection",
                            "symbols": [
                                "#main/Multiple_Alignment_Gene_Selection/Multiple_Alignment_Gene_Selection/Any Overlap",
                                "#main/Multiple_Alignment_Gene_Selection/Multiple_Alignment_Gene_Selection/Most Base Overlap"
                            ]
                        }
                    ],
                    "id": "#main/Multiple_Alignment_Gene_Selection"
                },
                {
                    "label": "Multiple Alignment Handling",
                    "doc": "Specify how to handle multiple alignments. By default reads with multiple alignments are discarded.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Multiple_Alignment_Handling/Multiple_Alignment_Handling",
                            "symbols": [
                                "#main/Multiple_Alignment_Handling/Multiple_Alignment_Handling/Discard",
                                "#main/Multiple_Alignment_Handling/Multiple_Alignment_Handling/Uniform Random"
                            ]
                        }
                    ],
                    "id": "#main/Multiple_Alignment_Handling"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "ATAC Predefined Peak Regions",
                    "doc": "An optional BED file containing pre-established chromatin accessibility peak regions for generating the ATAC cell-by-peak matrix. Only applies to ATAC assays.",
                    "id": "#main/Predefined_ATAC_Peaks"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "Reads",
                    "doc": "FASTQ files from libraries that may include WTA mRNA, Targeted mRNA, AbSeq, Sample Multiplexing, and related technologies",
                    "id": "#main/Reads"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "Reads-ATAC",
                    "doc": "FASTQ files from libraries generated using the ATAC assay protocol. Each lane of a library is expected to have 3 FASTQs - R1, R2 and I1/I2, where the index read contains the Cell Barcode and UMI sequence. Only applies to ATAC assays.",
                    "id": "#main/Reads_ATAC"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Reference Files Archive",
                    "id": "#main/Reference_Archive"
                },
                {
                    "label": "Run Name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "This is a name for output files, for example Experiment1_Metrics_Summary.csv. Default if left empty is to name run based on a library. Any non-alpha numeric characters will be changed to a hyphen.",
                    "id": "#main/Run_Name"
                },
                {
                    "label": "Sample Tags Version",
                    "doc": "The sample multiplexing kit version.  This option should only be set for a multiplexed experiment.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/Sample_Tags_Version/Sample_Tags_Version",
                            "symbols": [
                                "#main/Sample_Tags_Version/Sample_Tags_Version/human",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/hs",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/mouse",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/mm",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/flex",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/nuclei_includes_mrna",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/nuclei_atac_only",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/custom"
                            ]
                        }
                    ],
                    "id": "#main/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "Supplemental Reference",
                    "id": "#main/Supplemental_Reference"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "label": "Sample Tag Names",
                    "doc": "Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Should be alpha numeric, with + - and _ allowed. Any special characters: &, (), [], {}, <>, ?, | will be corrected to underscores. \n",
                    "id": "#main/Tag_Names"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#main/Target_analysis"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "Targeted Reference",
                    "id": "#main/Targeted_Reference"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "e-value threshold for J gene",
                    "doc": "The e-value threshold for J gene call by IgBlast, default is set as 0.001\n",
                    "id": "#main/VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "e-value threshold for V gene",
                    "doc": "The e-value threshold for V gene call by IgBlast, default is set as 0.001\n",
                    "id": "#main/VDJ_VGene_Evalue"
                },
                {
                    "label": "VDJ Species Version",
                    "doc": "The VDJ species and chain types.  This option should only be set for VDJ experiment.",
                    "type": [
                        "null",
                        {
                            "type": "enum",
                            "name": "#main/VDJ_Version/VDJ_Version",
                            "symbols": [
                                "#main/VDJ_Version/VDJ_Version/human",
                                "#main/VDJ_Version/VDJ_Version/hs",
                                "#main/VDJ_Version/VDJ_Version/mouse",
                                "#main/VDJ_Version/VDJ_Version/mm",
                                "#main/VDJ_Version/VDJ_Version/humanBCR",
                                "#main/VDJ_Version/VDJ_Version/humanTCR",
                                "#main/VDJ_Version/VDJ_Version/mouseBCR",
                                "#main/VDJ_Version/VDJ_Version/mouseTCR"
                            ]
                        }
                    ],
                    "id": "#main/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#main/Write_Filtered_Reads"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/MergeATAC/ATAC_out",
                    "id": "#main/ATAC"
                },
                {
                    "label": "BAM files and indices",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": [
                        "#main/MergeBAM_RNA/Bam",
                        "#main/MergeBAM_RNA/BamIndex",
                        "#main/MergeBAM_VDJ/Bam",
                        "#main/MergeBAM_VDJ/BamIndex",
                        "#main/MergeBAM_ATAC/Bam",
                        "#main/MergeBAM_ATAC/BamIndex"
                    ],
                    "linkMerge": "merge_flattened",
                    "pickValue": "all_non_null",
                    "id": "#main/Bam"
                },
                {
                    "label": "Bioproduct Statistics",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Bioproduct_Stats",
                    "id": "#main/Bioproduct_Stats"
                },
                {
                    "label": "Biotype Counts",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GenerateH5MU/biotypeCounts",
                    "id": "#main/BiotypeCounts"
                },
                {
                    "label": "Data Tables",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/GetDataTable/Data_Tables",
                    "id": "#main/Data_Tables"
                },
                {
                    "label": "Dimensionality Reduction Coordinates",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": [
                        "#main/GetDataTable/Dim_Reduction_Coord",
                        "#main/ATAC_Compile_Results/ATAC_Dim_Reduction_Coord"
                    ],
                    "pickValue": "all_non_null",
                    "linkMerge": "merge_flattened",
                    "id": "#main/Dim_Reduction_Coord"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/QualCLAlign_ATAC/Failed_Reads_CSVs",
                    "id": "#main/Failed_Reads_CSVs_ATAC"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/QualCLAlign_RNA/Failed_Reads_CSVs",
                    "id": "#main/Failed_Reads_CSVs_RNA"
                },
                {
                    "label": "Cellismo File",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GenerateH5MU/H5MU",
                    "id": "#main/H5MU"
                },
                {
                    "label": "Immune Cell Classification (Experimental)",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": [
                        "#main/GetDataTable/Cell_Type_Predictions",
                        "#main/ATAC_Compile_Results/Cell_Type_Predictions"
                    ],
                    "pickValue": "all_non_null",
                    "linkMerge": "merge_flattened",
                    "id": "#main/Immune_Cell_Classification(Experimental)"
                },
                {
                    "label": "Pipeline Logs",
                    "type": "Directory",
                    "outputSource": "#main/BundleLogs/logs_dir",
                    "id": "#main/Logs"
                },
                {
                    "label": "Metrics Summary",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/Metrics/Metrics_Summary",
                    "id": "#main/Metrics_Summary"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/MergeMultiplex/Multiplex_out",
                    "id": "#main/Multiplexing"
                },
                {
                    "label": "Pipeline Report HTML",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": [
                        "#main/GenerateH5MU/Pipeline_Report_HTML_Inj",
                        "#main/Metrics/Pipeline_Report_HTML"
                    ],
                    "pickValue": "first_non_null",
                    "id": "#main/Pipeline_Report_HTML_final"
                },
                {
                    "label": "Protein Aggregates (Experimental)",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GetDataTable/Protein_Aggregates_Experimental",
                    "id": "#main/Protein_Aggregates_Experimental"
                },
                {
                    "label": "Seurat RDS File",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/GenerateSeurat/SeuratRDS",
                    "id": "#main/Seurat"
                },
                {
                    "label": "vdjCellsDatatable",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjCellsDatatable",
                    "id": "#main/vdjCellsDatatable"
                },
                {
                    "label": "vdjCellsDatatableUncorrected",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjCellsDatatableUncorrected",
                    "id": "#main/vdjCellsDatatableUncorrected"
                },
                {
                    "label": "vdjDbecFilterImages",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjDbecFilterImages",
                    "id": "#main/vdjDbecFilterImages"
                },
                {
                    "label": "vdjDominantContigsAIRR",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjDominantContigsAIRR",
                    "id": "#main/vdjDominantContigsAIRR"
                },
                {
                    "label": "vdjMetricsCsv",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjMetricsCsv",
                    "id": "#main/vdjMetricsCsv"
                },
                {
                    "label": "vdjUnfilteredContigsAIRR",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#main/VDJ_Compile_Results/vdjUnfilteredContigsAIRR",
                    "id": "#main/vdjUnfilteredContigsAIRR"
                }
            ],
            "steps": [
                {
                    "run": "#ATAC_Cell_by_Peak.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_ATAC",
                            "id": "#main/ATAC_Cell_by_Peak/Assay"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Fragments",
                            "id": "#main/ATAC_Cell_by_Peak/Fragments"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/ATAC_Cell_by_Peak/GTF"
                        },
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "valueFrom": "$(self == null)",
                            "id": "#main/ATAC_Cell_by_Peak/Map_Gene_Name_ID"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Peaks",
                            "id": "#main/ATAC_Cell_by_Peak/Peaks"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/ATAC_Cell_by_Peak/Run_Metadata"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Transposase_Sites",
                            "id": "#main/ATAC_Cell_by_Peak/Transposase_Sites"
                        }
                    ],
                    "out": [
                        "#main/ATAC_Cell_by_Peak/Initial_Seurat_RDS",
                        "#main/ATAC_Cell_by_Peak/ATAC_Cell_Calling_Data",
                        "#main/ATAC_Cell_by_Peak/Total_Fragment_Metrics",
                        "#main/ATAC_Cell_by_Peak/Gene_Name_ID_Json",
                        "#main/ATAC_Cell_by_Peak/output"
                    ],
                    "when": "$(inputs.Assay == \"ATAC\")",
                    "id": "#main/ATAC_Cell_by_Peak"
                },
                {
                    "run": "#ATAC_Compile_Results.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_ATAC",
                            "id": "#main/ATAC_Compile_Results/Assay"
                        },
                        {
                            "source": "#main/GetDataTable/Biop_putative_data_table",
                            "id": "#main/ATAC_Compile_Results/Biop_putative_data_table"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/ATAC_Compile_Results/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Order_Subsampled",
                            "id": "#main/ATAC_Compile_Results/Cell_Order_Subsampled"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Fragments",
                            "id": "#main/ATAC_Compile_Results/Fragments"
                        },
                        {
                            "source": [
                                "#main/AlignmentAnalysis/Gene_Name_ID_Json",
                                "#main/ATAC_Cell_by_Peak/Gene_Name_ID_Json"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/ATAC_Compile_Results/Gene_Name_ID_Json"
                        },
                        {
                            "source": "#main/ATAC_Cell_by_Peak/Initial_Seurat_RDS",
                            "id": "#main/ATAC_Compile_Results/Initial_Seurat_RDS"
                        },
                        {
                            "source": "#main/GetDataTable/Metrics_tar",
                            "id": "#main/ATAC_Compile_Results/Input_Metrics_tar"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Reference_Genome_Contig_Lengths",
                            "id": "#main/ATAC_Compile_Results/Reference_Genome_Contig_Lengths"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Reference_Genome_Size",
                            "id": "#main/ATAC_Compile_Results/Reference_Genome_Size"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/ATAC_Compile_Results/Run_Metadata"
                        },
                        {
                            "source": "#main/ATAC_Cell_by_Peak/Total_Fragment_Metrics",
                            "id": "#main/ATAC_Compile_Results/Total_Fragment_Metrics"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/UnifiedMetrics",
                            "id": "#main/ATAC_Compile_Results/Unified_Metrics"
                        }
                    ],
                    "out": [
                        "#main/ATAC_Compile_Results/Metrics_tar",
                        "#main/ATAC_Compile_Results/ATAC_Data_Tables",
                        "#main/ATAC_Compile_Results/ATAC_Dim_Reduction_Coord",
                        "#main/ATAC_Compile_Results/ATAC_Seurat_RDS",
                        "#main/ATAC_Compile_Results/Cell_Type_Predictions",
                        "#main/ATAC_Compile_Results/ATAC_TSS_Enrichment",
                        "#main/ATAC_Compile_Results/ATAC_Nucleosome_Signal",
                        "#main/ATAC_Compile_Results/ATAC_Transposase_Site_Counts",
                        "#main/ATAC_Compile_Results/output"
                    ],
                    "when": "$(inputs.Assay == \"ATAC\")",
                    "id": "#main/ATAC_Compile_Results"
                },
                {
                    "run": "#AddtoBam.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_ATAC",
                            "id": "#main/AddtoBam_ATAC/Assay"
                        },
                        {
                            "source": [
                                "#main/QualCLAlign_ATAC/BAMFiles"
                            ],
                            "default": [
                                "does_not_exist"
                            ],
                            "id": "#main/AddtoBam_ATAC/Bam"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/AddtoBam_ATAC/Cell_Order"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/AddtoBam_ATAC/Generate_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AddtoBam_ATAC/Run_Metadata"
                        },
                        {
                            "source": "#main/SampleTag/SampleTag_Calls",
                            "id": "#main/AddtoBam_ATAC/SampleTag_Calls"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true && inputs.Assay == \"ATAC\")",
                    "out": [
                        "#main/AddtoBam_ATAC/Annotated_Bam",
                        "#main/AddtoBam_ATAC/output"
                    ],
                    "scatter": [
                        "#main/AddtoBam_ATAC/Bam"
                    ],
                    "id": "#main/AddtoBam_ATAC"
                },
                {
                    "run": "#AddtoBam.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/AddtoBam_RNA/Assay"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Annotated_Bam_Files",
                            "default": [
                                "does_not_exist"
                            ],
                            "id": "#main/AddtoBam_RNA/Bam"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/AddtoBam_RNA/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/AddtoBam_RNA/Cell_Order"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Corrected_Mols_List",
                            "id": "#main/AddtoBam_RNA/Corrected_Mols"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/AddtoBam_RNA/Generate_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AddtoBam_RNA/Run_Metadata"
                        },
                        {
                            "source": "#main/SampleTag/SampleTag_Calls",
                            "id": "#main/AddtoBam_RNA/SampleTag_Calls"
                        },
                        {
                            "source": "#main/CheckReference/Target_Gene_Mapping",
                            "id": "#main/AddtoBam_RNA/Target_Gene_Mapping"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true && inputs.Cell_Calling_Data != 6 && (inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\"))",
                    "out": [
                        "#main/AddtoBam_RNA/Annotated_Bam",
                        "#main/AddtoBam_RNA/output"
                    ],
                    "scatter": [
                        "#main/AddtoBam_RNA/Bam"
                    ],
                    "id": "#main/AddtoBam_RNA"
                },
                {
                    "run": "#AddtoBam.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/AddtoBam_VDJ/Assay"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Annotated_Bam_Files",
                            "default": [
                                "does_not_exist"
                            ],
                            "id": "#main/AddtoBam_VDJ/Bam"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/AddtoBam_VDJ/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/VDJ_Compile_Results/Cell_Order",
                            "id": "#main/AddtoBam_VDJ/Cell_Order"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Corrected_Mols_List",
                            "id": "#main/AddtoBam_VDJ/Corrected_Mols"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/AddtoBam_VDJ/Generate_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AddtoBam_VDJ/Run_Metadata"
                        },
                        {
                            "source": "#main/SampleTag/SampleTag_Calls",
                            "id": "#main/AddtoBam_VDJ/SampleTag_Calls"
                        },
                        {
                            "source": "#main/CheckReference/Target_Gene_Mapping",
                            "id": "#main/AddtoBam_VDJ/Target_Gene_Mapping"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true && inputs.Cell_Calling_Data == 6 && (inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\"))",
                    "out": [
                        "#main/AddtoBam_VDJ/Annotated_Bam",
                        "#main/AddtoBam_VDJ/output"
                    ],
                    "scatter": [
                        "#main/AddtoBam_VDJ/Bam"
                    ],
                    "id": "#main/AddtoBam_VDJ"
                },
                {
                    "run": "#AlignmentAnalysis.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/AlignmentAnalysis/Assay"
                        },
                        {
                            "source": "#main/Misc_Settings/Exclude_Intronic_Reads",
                            "id": "#main/AlignmentAnalysis/Exclude_Intronic_Reads"
                        },
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AlignmentAnalysis/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/AlignmentAnalysis/GTF"
                        },
                        {
                            "source": "#main/Maximum_Threads",
                            "id": "#main/AlignmentAnalysis/Maximum_Threads"
                        },
                        {
                            "source": "#main/Misc_Settings/Multiple_Alignment_Gene_Selection",
                            "id": "#main/AlignmentAnalysis/Multiple_Alignment_Gene_Selection"
                        },
                        {
                            "source": "#main/Misc_Settings/Multiple_Alignment_Handling",
                            "id": "#main/AlignmentAnalysis/Multiple_Alignment_Handling"
                        },
                        {
                            "source": "#main/QualCLAlign_RNA/BAMFiles",
                            "id": "#main/AlignmentAnalysis/R2_BAM"
                        },
                        {
                            "source": "#main/QualCLAlign_RNA/QualCLAlignMetrics",
                            "id": "#main/AlignmentAnalysis/ReadQualityMetrics"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AlignmentAnalysis/Run_Metadata"
                        },
                        {
                            "source": "#main/CheckReference/Transcript_Length",
                            "id": "#main/AlignmentAnalysis/Transcript_Length"
                        }
                    ],
                    "out": [
                        "#main/AlignmentAnalysis/Seq_Metrics",
                        "#main/AlignmentAnalysis/Annotated_Bam_Files",
                        "#main/AlignmentAnalysis/Sorted_Valid_Reads_CSV",
                        "#main/AlignmentAnalysis/Gene_Name_ID_Json",
                        "#main/AlignmentAnalysis/num_valid_ig_reads",
                        "#main/AlignmentAnalysis/num_valid_tcr_reads",
                        "#main/AlignmentAnalysis/validIgReads",
                        "#main/AlignmentAnalysis/validTcrReads",
                        "#main/AlignmentAnalysis/num_cell_estimate",
                        "#main/AlignmentAnalysis/num_bioproducts"
                    ],
                    "when": "$(inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\")",
                    "id": "#main/AlignmentAnalysis"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 32000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#AnnotateMolecules.cwl",
                    "in": [
                        {
                            "source": "#main/Internal_Settings/AbSeq_UMI",
                            "id": "#main/AnnotateMolecules/AbSeq_UMI"
                        },
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/AnnotateMolecules/Assay"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/AnnotateMolecules/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateMolecules/Run_Metadata"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/AnnotateMolecules/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Sorted_Valid_Reads_CSV",
                            "default": [
                                "does_not_exist"
                            ],
                            "id": "#main/AnnotateMolecules/Valids"
                        }
                    ],
                    "out": [
                        "#main/AnnotateMolecules/Bioproduct_Stats_List",
                        "#main/AnnotateMolecules/Cell_Biop_Summary_List",
                        "#main/AnnotateMolecules/Corrected_Mols_List",
                        "#main/AnnotateMolecules/Total_Molecules",
                        "#main/AnnotateMolecules/Num_Molecules_Corrected",
                        "#main/AnnotateMolecules/output"
                    ],
                    "scatter": [
                        "#main/AnnotateMolecules/Valids"
                    ],
                    "when": "$((inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\") && (inputs.Cell_Calling_Data != 6 || inputs.Sample_Tags_Version != null))",
                    "id": "#main/AnnotateMolecules"
                },
                {
                    "run": "#Assay_Settings.cwl",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "valueFrom": "${ if (self && self.length > 0){ return true; } else{ return false; } }",
                            "id": "#main/Assay_Settings/AbSeq_Reference_Present"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/Assay_Settings/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/Reads_ATAC",
                            "valueFrom": "${ if (self && self.length > 0) { return true; } else { return false; } }",
                            "id": "#main/Assay_Settings/Reads_ATAC_Present"
                        },
                        {
                            "source": "#main/Reads",
                            "valueFrom": "${ if (self && self.length > 0) { return true; } else { return false; } }",
                            "id": "#main/Assay_Settings/Reads_RNA_Present"
                        },
                        {
                            "source": "#main/Reference_Archive",
                            "valueFrom": "${ if (self){ return true; } else{ return false; } }",
                            "id": "#main/Assay_Settings/Reference_Archive_Present"
                        },
                        {
                            "source": "#main/Targeted_Reference",
                            "valueFrom": "${ if (self && self.length > 0){ return true; } else{ return false; } }",
                            "id": "#main/Assay_Settings/Targeted_Reference_Present"
                        }
                    ],
                    "out": [
                        "#main/Assay_Settings/Assay_RNA",
                        "#main/Assay_Settings/Assay_ATAC"
                    ],
                    "id": "#main/Assay_Settings"
                },
                {
                    "run": "#ATAC_Motif_Analysis.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_ATAC",
                            "id": "#main/AtacMotifAnalysis/Assay"
                        },
                        {
                            "source": "#main/ATAC_Compile_Results/ATAC_Data_Tables",
                            "id": "#main/AtacMotifAnalysis/Atac_Datatables"
                        },
                        {
                            "source": "#main/Maximum_Threads",
                            "id": "#main/AtacMotifAnalysis/Maximum_Threads"
                        },
                        {
                            "source": "#main/CheckReference/Reference_FASTA",
                            "id": "#main/AtacMotifAnalysis/Reference_FASTA"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AtacMotifAnalysis/Run_Metadata"
                        },
                        {
                            "source": "#main/CheckReference/Transcription_Factor_Motif_PFM",
                            "id": "#main/AtacMotifAnalysis/Transcription_Factor_Motif_PFM"
                        }
                    ],
                    "out": [
                        "#main/AtacMotifAnalysis/ATAC_Motif_Data_Tables",
                        "#main/AtacMotifAnalysis/ATAC_Peak_Seq",
                        "#main/AtacMotifAnalysis/ATAC_Background_Peaks",
                        "#main/AtacMotifAnalysis/ATAC_Motif_Positions",
                        "#main/AtacMotifAnalysis/output"
                    ],
                    "when": "$(inputs.Assay == \"ATAC\")",
                    "id": "#main/AtacMotifAnalysis"
                },
                {
                    "label": "Bam Settings",
                    "run": "#BamSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Generate_Bam",
                            "id": "#main/Bam_Settings/_Generate_Bam"
                        }
                    ],
                    "out": [
                        "#main/Bam_Settings/Generate_Bam"
                    ],
                    "id": "#main/Bam_Settings"
                },
                {
                    "run": "#BundleLogs.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/CheckReference/output",
                                "#main/GetDataTable/output",
                                "#main/ATAC_Cell_by_Peak/output",
                                "#main/ATAC_Compile_Results/output",
                                "#main/Metrics/output",
                                "#main/AddtoBam_RNA/output",
                                "#main/AddtoBam_VDJ/output",
                                "#main/AddtoBam_ATAC/output",
                                "#main/AnnotateMolecules/output",
                                "#main/MergeBAM_RNA/log",
                                "#main/MergeBAM_VDJ/log",
                                "#main/MergeBAM_ATAC/log",
                                "#main/GenerateH5MU/output",
                                "#main/GenerateSeurat/output",
                                "#main/Peak_Annotation/output",
                                "#main/CalcGeneBodyCoverage/output"
                            ],
                            "pickValue": "all_non_null",
                            "linkMerge": "merge_flattened",
                            "id": "#main/BundleLogs/log_files"
                        }
                    ],
                    "out": [
                        "#main/BundleLogs/logs_dir"
                    ],
                    "id": "#main/BundleLogs"
                },
                {
                    "run": "#CalcGeneBodyCoverage.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/CalcGeneBodyCoverage/Assay"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Annotated_Bam_Files",
                            "id": "#main/CalcGeneBodyCoverage/Bam"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/CalcGeneBodyCoverage/Gtf"
                        },
                        {
                            "source": "#main/Maximum_Threads",
                            "id": "#main/CalcGeneBodyCoverage/Maximum_Threads"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Base_Name",
                            "id": "#main/CalcGeneBodyCoverage/Run_Base_Name"
                        }
                    ],
                    "when": "$(inputs.Assay == \"WTA\")",
                    "out": [
                        "#main/CalcGeneBodyCoverage/Gene_Body_Coverage_Json",
                        "#main/CalcGeneBodyCoverage/output"
                    ],
                    "id": "#main/CalcGeneBodyCoverage"
                },
                {
                    "run": "#CalcSeqSaturation.cwl",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "valueFrom": "${ if (self && self.length > 0){ return true; } else{ return false; } }",
                            "id": "#main/CalcSeqSaturation/AbSeq_Reference_Present"
                        },
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/CalcSeqSaturation/Assay"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/CalcSeqSaturation/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Corrected_Mols_List",
                            "id": "#main/CalcSeqSaturation/Corrected_Mols"
                        },
                        {
                            "source": "#main/GetDataTable/Metrics_tar",
                            "id": "#main/CalcSeqSaturation/Metrics_Tar"
                        },
                        {
                            "source": "#main/SumMoleculesCorrected/Num_Mols_Corrected",
                            "id": "#main/CalcSeqSaturation/Num_Mols_Corrected"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/CalcSeqSaturation/Putative_Cells"
                        },
                        {
                            "source": "#main/Reference_Archive",
                            "valueFrom": "${ if (self){ return true; } else{ return false; } }",
                            "id": "#main/CalcSeqSaturation/Reference_Archive_Present"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Base_Name",
                            "id": "#main/CalcSeqSaturation/Run_Base_Name"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/CalcSeqSaturation/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Targeted_Reference",
                            "valueFrom": "${ if (self && self.length > 0){ return true; } else{ return false; } }",
                            "id": "#main/CalcSeqSaturation/Targeted_Reference_Present"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Sorted_Valid_Reads_CSV",
                            "id": "#main/CalcSeqSaturation/Valid_Reads_CSV"
                        }
                    ],
                    "out": [
                        "#main/CalcSeqSaturation/Seq_Saturation_Plot_Json"
                    ],
                    "when": "$((inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\") && inputs.Cell_Calling_Data != 6 && inputs.Sample_Tags_Version != \"nuclei_atac_only\")",
                    "id": "#main/CalcSeqSaturation"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 10000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#CheckReference.cwl",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "id": "#main/CheckReference/AbSeq_Reference"
                        },
                        {
                            "source": [
                                "#main/Assay_Settings/Assay_RNA",
                                "#main/Assay_Settings/Assay_ATAC"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/CheckReference/Assay"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/CheckReference/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/Predefined_ATAC_Peaks",
                            "id": "#main/CheckReference/Predefined_ATAC_Peaks"
                        },
                        {
                            "source": "#main/Reference_Archive",
                            "id": "#main/CheckReference/Reference_Archive"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/CheckReference/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Supplemental_Reference",
                            "id": "#main/CheckReference/Supplemental_Reference"
                        },
                        {
                            "source": "#main/Targeted_Reference",
                            "id": "#main/CheckReference/Targeted_Reference"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/CheckReference/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#main/CheckReference/Index",
                        "#main/CheckReference/Extra_Seqs",
                        "#main/CheckReference/Full_Genes",
                        "#main/CheckReference/output",
                        "#main/CheckReference/Transcript_Length",
                        "#main/CheckReference/GTF",
                        "#main/CheckReference/Reference_FASTA",
                        "#main/CheckReference/Target_Gene_Mapping",
                        "#main/CheckReference/Transcription_Factor_Motif_PFM",
                        "#main/CheckReference/Checked_Predefined_Peaks",
                        "#main/CheckReference/Checked_Predefined_Peaks_Index"
                    ],
                    "id": "#main/CheckReference"
                },
                {
                    "run": "#GenerateH5MU.cwl",
                    "in": [
                        {
                            "source": "#main/AtacMotifAnalysis/ATAC_Background_Peaks",
                            "id": "#main/GenerateH5MU/Atac_Background_Peaks"
                        },
                        {
                            "source": [
                                "#main/ATAC_Compile_Results/ATAC_Data_Tables",
                                "#main/AtacMotifAnalysis/ATAC_Motif_Data_Tables"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/GenerateH5MU/Atac_Datatables"
                        },
                        {
                            "source": "#main/Metrics/Metrics_ATAC",
                            "id": "#main/GenerateH5MU/Atac_Metrics"
                        },
                        {
                            "source": "#main/AtacMotifAnalysis/ATAC_Motif_Positions",
                            "id": "#main/GenerateH5MU/Atac_Motif_Positions"
                        },
                        {
                            "source": "#main/ATAC_Compile_Results/ATAC_Nucleosome_Signal",
                            "id": "#main/GenerateH5MU/Atac_Nucleosome_Signal"
                        },
                        {
                            "source": "#main/AtacMotifAnalysis/ATAC_Peak_Seq",
                            "id": "#main/GenerateH5MU/Atac_Peak_Seq"
                        },
                        {
                            "source": "#main/ATAC_Compile_Results/ATAC_TSS_Enrichment",
                            "id": "#main/GenerateH5MU/Atac_TSS_Enrichment"
                        },
                        {
                            "source": "#main/ATAC_Compile_Results/ATAC_Transposase_Site_Counts",
                            "id": "#main/GenerateH5MU/Atac_Transposase_Site_Counts"
                        },
                        {
                            "source": "#main/GetDataTable/Bioproduct_Stats",
                            "id": "#main/GenerateH5MU/Bioproduct_Stats"
                        },
                        {
                            "source": [
                                "#main/GetDataTable/Cell_Type_Predictions",
                                "#main/ATAC_Compile_Results/Cell_Type_Predictions"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/GenerateH5MU/Cell_Type_Experimental"
                        },
                        {
                            "source": "#main/GetDataTable/Data_Tables",
                            "id": "#main/GenerateH5MU/Data_Tables"
                        },
                        {
                            "source": [
                                "#main/GetDataTable/Dim_Reduction_Coord",
                                "#main/ATAC_Compile_Results/ATAC_Dim_Reduction_Coord"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/GenerateH5MU/Dim_Reduction_Coord"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/GenerateH5MU/GTF"
                        },
                        {
                            "source": "#main/Metrics/Metrics_Summary",
                            "id": "#main/GenerateH5MU/Metrics_Summary"
                        },
                        {
                            "source": "#main/Peak_Annotation/Peak_Annotation_TSV",
                            "id": "#main/GenerateH5MU/Peak_Annotation"
                        },
                        {
                            "source": "#main/Metrics/Pipeline_Report_HTML",
                            "id": "#main/GenerateH5MU/Pipeline_Report_Html"
                        },
                        {
                            "source": "#main/GetDataTable/Protein_Aggregates_Experimental",
                            "id": "#main/GenerateH5MU/Protein_Aggregates_Experimental"
                        },
                        {
                            "source": "#main/GetDataTable/Putative_Cells_Origin",
                            "id": "#main/GenerateH5MU/Putative_Cells_Origin"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Reference_Genome_Contig_Lengths",
                            "id": "#main/GenerateH5MU/Reference_Genome_Contig_Lengths"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/GenerateH5MU/Run_Metadata"
                        },
                        {
                            "source": "#main/SampleTag/SampleTag_CSVs",
                            "id": "#main/GenerateH5MU/SampleTag_CSVs"
                        },
                        {
                            "source": "#main/SampleTag/SampleTag_Calls",
                            "id": "#main/GenerateH5MU/SampleTag_Calls"
                        },
                        {
                            "source": "#main/CheckReference/Transcription_Factor_Motif_PFM",
                            "id": "#main/GenerateH5MU/Transcription_Factor_Motif_PFM"
                        },
                        {
                            "source": "#main/VDJ_Compile_Results/vdjAnnData",
                            "id": "#main/GenerateH5MU/VDJ_AnnData"
                        },
                        {
                            "source": "#main/VDJ_Compile_Results/vdjCellsDatatable",
                            "id": "#main/GenerateH5MU/VDJ_Per_Cell"
                        }
                    ],
                    "out": [
                        "#main/GenerateH5MU/H5MU",
                        "#main/GenerateH5MU/output",
                        "#main/GenerateH5MU/biotypeCounts",
                        "#main/GenerateH5MU/Pipeline_Report_HTML_Inj"
                    ],
                    "id": "#main/GenerateH5MU"
                },
                {
                    "run": "#GenerateSeurat.cwl",
                    "in": [
                        {
                            "source": "#main/ATAC_Compile_Results/ATAC_Seurat_RDS",
                            "id": "#main/GenerateSeurat/ATAC_Seurat"
                        },
                        {
                            "source": "#main/AtacMotifAnalysis/ATAC_Background_Peaks",
                            "id": "#main/GenerateSeurat/Atac_Background_Peaks"
                        },
                        {
                            "source": "#main/ATAC_Compile_Results/ATAC_Data_Tables",
                            "id": "#main/GenerateSeurat/Atac_Datatables"
                        },
                        {
                            "source": "#main/Metrics/Metrics_ATAC",
                            "id": "#main/GenerateSeurat/Atac_Metrics"
                        },
                        {
                            "source": "#main/AtacMotifAnalysis/ATAC_Motif_Data_Tables",
                            "id": "#main/GenerateSeurat/Atac_Motif_Datatables"
                        },
                        {
                            "source": "#main/AtacMotifAnalysis/ATAC_Motif_Positions",
                            "id": "#main/GenerateSeurat/Atac_Motif_Positions"
                        },
                        {
                            "source": "#main/AtacMotifAnalysis/ATAC_Peak_Seq",
                            "id": "#main/GenerateSeurat/Atac_Peak_Seq"
                        },
                        {
                            "source": "#main/GetDataTable/Bioproduct_Stats",
                            "id": "#main/GenerateSeurat/Bioproduct_Stats"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/GenerateSeurat/Cell_Calling_Data"
                        },
                        {
                            "source": [
                                "#main/GetDataTable/Cell_Type_Predictions",
                                "#main/ATAC_Compile_Results/Cell_Type_Predictions"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/GenerateSeurat/Cell_Type_Experimental"
                        },
                        {
                            "source": "#main/GetDataTable/Data_Tables",
                            "id": "#main/GenerateSeurat/Data_Tables"
                        },
                        {
                            "source": [
                                "#main/GetDataTable/Dim_Reduction_Coord",
                                "#main/ATAC_Compile_Results/ATAC_Dim_Reduction_Coord"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/GenerateSeurat/Dim_Reduction_Coord"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/GenerateSeurat/GTF"
                        },
                        {
                            "source": "#main/Metrics/Metrics_Summary",
                            "id": "#main/GenerateSeurat/Metrics_Summary"
                        },
                        {
                            "source": "#main/Peak_Annotation/Peak_Annotation_TSV",
                            "id": "#main/GenerateSeurat/Peak_Annotation"
                        },
                        {
                            "source": "#main/GetDataTable/Protein_Aggregates_Experimental",
                            "id": "#main/GenerateSeurat/Protein_Aggregates_Experimental"
                        },
                        {
                            "source": "#main/GetDataTable/Putative_Cells_Origin",
                            "id": "#main/GenerateSeurat/Putative_Cells_Origin"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Reference_Genome_Contig_Lengths",
                            "id": "#main/GenerateSeurat/Reference_Genome_Contig_Lengths"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/GenerateSeurat/Run_Metadata"
                        },
                        {
                            "source": "#main/SampleTag/SampleTag_CSVs",
                            "id": "#main/GenerateSeurat/SampleTag_CSVs"
                        },
                        {
                            "source": "#main/SampleTag/SampleTag_Calls",
                            "id": "#main/GenerateSeurat/SampleTag_Calls"
                        },
                        {
                            "source": "#main/CheckReference/Transcription_Factor_Motif_PFM",
                            "id": "#main/GenerateSeurat/Transcription_Factor_Motif_PFM"
                        },
                        {
                            "source": "#main/VDJ_Compile_Results/vdjCellsDatatable",
                            "id": "#main/GenerateSeurat/VDJ_Per_Cell"
                        }
                    ],
                    "out": [
                        "#main/GenerateSeurat/SeuratRDS",
                        "#main/GenerateSeurat/output"
                    ],
                    "when": "$(inputs.Cell_Calling_Data != 6)",
                    "id": "#main/GenerateSeurat"
                },
                {
                    "run": "#GetDataTable.cwl",
                    "in": [
                        {
                            "source": "#main/ATAC_Cell_by_Peak/ATAC_Cell_Calling_Data",
                            "id": "#main/GetDataTable/ATAC_Cell_Calling_Input"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Bioproduct_Stats_List",
                            "pickValue": "all_non_null",
                            "id": "#main/GetDataTable/Bioproduct_Stats_List"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Cell_Biop_Summary_List",
                            "pickValue": "all_non_null",
                            "id": "#main/GetDataTable/Cell_Biop_Summary_List"
                        },
                        {
                            "source": "#main/CheckReference/Full_Genes",
                            "id": "#main/GetDataTable/Full_Genes"
                        },
                        {
                            "source": [
                                "#main/AlignmentAnalysis/Gene_Name_ID_Json",
                                "#main/ATAC_Cell_by_Peak/Gene_Name_ID_Json"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/GetDataTable/Gene_Name_ID_Json"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/GetDataTable/Run_Metadata"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/Seq_Metrics",
                            "id": "#main/GetDataTable/Seq_Metrics"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Total_Molecules",
                            "pickValue": "all_non_null",
                            "id": "#main/GetDataTable/Total_Molecules"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_bioproducts",
                            "id": "#main/GetDataTable/num_bioproducts"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_cell_estimate",
                            "id": "#main/GetDataTable/num_cell_estimate"
                        }
                    ],
                    "out": [
                        "#main/GetDataTable/Metrics_tar",
                        "#main/GetDataTable/Bioproduct_Stats",
                        "#main/GetDataTable/Cell_Order",
                        "#main/GetDataTable/Cell_Order_Subsampled",
                        "#main/GetDataTable/Cell_Type_Predictions",
                        "#main/GetDataTable/Data_Tables",
                        "#main/GetDataTable/Dim_Reduction_Coord",
                        "#main/GetDataTable/output",
                        "#main/GetDataTable/Protein_Aggregates_Experimental",
                        "#main/GetDataTable/Putative_Cells_Origin",
                        "#main/GetDataTable/Biop_putative_data_table"
                    ],
                    "id": "#main/GetDataTable"
                },
                {
                    "label": "Get Machine Resources",
                    "run": "#GetMachineResources.cwl",
                    "in": [],
                    "out": [
                        "#main/GetMachineResources/Total_CPUs_Avail"
                    ],
                    "id": "#main/GetMachineResources"
                },
                {
                    "label": "Internal Settings",
                    "run": "#InternalSettings.cwl",
                    "in": [
                        {
                            "source": "#main/AbSeq_UMI",
                            "id": "#main/Internal_Settings/_AbSeq_UMI"
                        },
                        {
                            "source": "#main/Target_analysis",
                            "id": "#main/Internal_Settings/_Target_analysis"
                        },
                        {
                            "source": "#main/VDJ_JGene_Evalue",
                            "id": "#main/Internal_Settings/_VDJ_JGene_Evalue"
                        },
                        {
                            "source": "#main/VDJ_VGene_Evalue",
                            "id": "#main/Internal_Settings/_VDJ_VGene_Evalue"
                        }
                    ],
                    "out": [
                        "#main/Internal_Settings/AbSeq_UMI",
                        "#main/Internal_Settings/Target_analysis",
                        "#main/Internal_Settings/VDJ_VGene_Evalue",
                        "#main/Internal_Settings/VDJ_JGene_Evalue"
                    ],
                    "id": "#main/Internal_Settings"
                },
                {
                    "run": {
                        "cwlVersion": "v1.2",
                        "class": "ExpressionTool",
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "inputs": [
                            {
                                "type": {
                                    "items": [
                                        "null",
                                        "File"
                                    ],
                                    "type": "array"
                                },
                                "id": "#main/MergeATAC/run/ATAC_Files"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    {
                                        "type": "array",
                                        "items": "File"
                                    }
                                ],
                                "id": "#main/MergeATAC/run/ATAC_out"
                            }
                        ],
                        "expression": "${\n  var fp_array = [];\n  for (var i = 0; i < inputs.ATAC_Files.length; i++) {\n    var fp = inputs.ATAC_Files[i];\n    if (fp != null) {\n      fp_array.push(fp);\n    }\n  }\n  return({\"ATAC_out\": fp_array});\n}"
                    },
                    "in": [
                        {
                            "source": [
                                "#main/QualCLAlign_ATAC/Fragments",
                                "#main/QualCLAlign_ATAC/Transposase_Sites",
                                "#main/QualCLAlign_ATAC/Peaks",
                                "#main/Peak_Annotation/Peak_Annotation_TSV",
                                "#main/Metrics/Metrics_ATAC",
                                "#main/ATAC_Compile_Results/ATAC_Data_Tables",
                                "#main/AtacMotifAnalysis/ATAC_Motif_Data_Tables",
                                "#main/AtacMotifAnalysis/ATAC_Motif_Positions"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/MergeATAC/ATAC_Files"
                        },
                        {
                            "source": "#main/Assay_Settings/Assay_ATAC",
                            "id": "#main/MergeATAC/Assay"
                        }
                    ],
                    "out": [
                        "#main/MergeATAC/ATAC_out"
                    ],
                    "when": "$(inputs.Assay == \"ATAC\")",
                    "id": "#main/MergeATAC"
                },
                {
                    "run": "#MergeBAM.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_ATAC",
                            "id": "#main/MergeBAM_ATAC/Assay"
                        },
                        {
                            "source": "#main/AddtoBam_ATAC/Annotated_Bam",
                            "id": "#main/MergeBAM_ATAC/BamFiles"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/MergeBAM_ATAC/Generate_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/MergeBAM_ATAC/Run_Metadata"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true && inputs.Assay == \"ATAC\")",
                    "out": [
                        "#main/MergeBAM_ATAC/Bam",
                        "#main/MergeBAM_ATAC/BamIndex",
                        "#main/MergeBAM_ATAC/log"
                    ],
                    "id": "#main/MergeBAM_ATAC"
                },
                {
                    "run": "#MergeBAM.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/MergeBAM_RNA/Assay"
                        },
                        {
                            "source": "#main/AddtoBam_RNA/Annotated_Bam",
                            "id": "#main/MergeBAM_RNA/BamFiles"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/MergeBAM_RNA/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/MergeBAM_RNA/Generate_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/MergeBAM_RNA/Run_Metadata"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true && inputs.Cell_Calling_Data != 6 && (inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\"))",
                    "out": [
                        "#main/MergeBAM_RNA/Bam",
                        "#main/MergeBAM_RNA/BamIndex",
                        "#main/MergeBAM_RNA/log"
                    ],
                    "id": "#main/MergeBAM_RNA"
                },
                {
                    "run": "#MergeBAM.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/MergeBAM_VDJ/Assay"
                        },
                        {
                            "source": "#main/AddtoBam_VDJ/Annotated_Bam",
                            "id": "#main/MergeBAM_VDJ/BamFiles"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/MergeBAM_VDJ/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/MergeBAM_VDJ/Generate_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/MergeBAM_VDJ/Run_Metadata"
                        }
                    ],
                    "when": "$(inputs.Generate_Bam == true && inputs.Cell_Calling_Data == 6 && (inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\"))",
                    "out": [
                        "#main/MergeBAM_VDJ/Bam",
                        "#main/MergeBAM_VDJ/BamIndex",
                        "#main/MergeBAM_VDJ/log"
                    ],
                    "id": "#main/MergeBAM_VDJ"
                },
                {
                    "run": {
                        "cwlVersion": "v1.2",
                        "class": "ExpressionTool",
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "inputs": [
                            {
                                "type": {
                                    "items": [
                                        "null",
                                        "File"
                                    ],
                                    "type": "array"
                                },
                                "id": "#main/MergeMultiplex/run/SampleTag_Files"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    {
                                        "type": "array",
                                        "items": "File"
                                    }
                                ],
                                "id": "#main/MergeMultiplex/run/Multiplex_out"
                            }
                        ],
                        "expression": "${\n  var fp_array = [];\n  for (var i = 0; i < inputs.SampleTag_Files.length; i++) {\n    var fp = inputs.SampleTag_Files[i];\n    if (fp != null) {\n      fp_array.push(fp);\n    }\n  }\n  return({\"Multiplex_out\": fp_array});\n}"
                    },
                    "in": [
                        {
                            "source": [
                                "#main/SampleTag/SampleTag_CSVs",
                                "#main/SampleTag/SampleTag_perTagZips"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/MergeMultiplex/SampleTag_Files"
                        }
                    ],
                    "out": [
                        "#main/MergeMultiplex/Multiplex_out"
                    ],
                    "id": "#main/MergeMultiplex"
                },
                {
                    "run": "#Metadata.cwl",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "id": "#main/Metadata_Settings/AbSeq_Reference"
                        },
                        {
                            "source": [
                                "#main/Assay_Settings/Assay_RNA",
                                "#main/Assay_Settings/Assay_ATAC"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/Metadata_Settings/Assay"
                        },
                        {
                            "source": [
                                "#main/QualCLAlign_RNA/Bead_Version",
                                "#main/QualCLAlign_ATAC/Bead_Version"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/Metadata_Settings/Bead_Version"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_ATAC_Algorithm",
                            "id": "#main/Metadata_Settings/Cell_Calling_ATAC_Algorithm"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Bioproduct_Algorithm",
                            "id": "#main/Metadata_Settings/Cell_Calling_Bioproduct_Algorithm"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/Metadata_Settings/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/Misc_Settings/Custom_STAR_Params",
                            "id": "#main/Metadata_Settings/Custom_STAR_Params"
                        },
                        {
                            "source": "#main/Misc_Settings/Custom_bwa_mem2_Params",
                            "id": "#main/Metadata_Settings/Custom_bwa_mem2_Params"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                            "id": "#main/Metadata_Settings/Exact_Cell_Count"
                        },
                        {
                            "source": "#main/Misc_Settings/Exclude_Intronic_Reads",
                            "id": "#main/Metadata_Settings/Exclude_Intronic_Reads"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Expected_Cell_Count",
                            "id": "#main/Metadata_Settings/Expected_Cell_Count"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/Metadata_Settings/Generate_Bam"
                        },
                        {
                            "source": "#main/QualCLAlign_RNA/Libraries",
                            "pickValue": "all_non_null",
                            "id": "#main/Metadata_Settings/Libraries"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Libraries",
                            "pickValue": "all_non_null",
                            "id": "#main/Metadata_Settings/Libraries_ATAC"
                        },
                        {
                            "source": "#main/Misc_Settings/Long_Reads",
                            "id": "#main/Metadata_Settings/Long_Reads"
                        },
                        {
                            "valueFrom": "BD Rhapsody Sequence Analysis Pipeline",
                            "id": "#main/Metadata_Settings/Pipeline_Name"
                        },
                        {
                            "source": "#main/Version/version",
                            "id": "#main/Metadata_Settings/Pipeline_Version"
                        },
                        {
                            "source": "#main/Predefined_ATAC_Peaks",
                            "id": "#main/Metadata_Settings/Predefined_ATAC_Peaks"
                        },
                        {
                            "source": "#main/QualCLAlign_RNA/ReadsList",
                            "pickValue": "all_non_null",
                            "id": "#main/Metadata_Settings/Reads"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/ReadsList",
                            "pickValue": "all_non_null",
                            "id": "#main/Metadata_Settings/Reads_ATAC"
                        },
                        {
                            "source": "#main/Reference_Archive",
                            "valueFrom": "${\n  return self ? self.basename : null;\n}\n",
                            "id": "#main/Metadata_Settings/Reference_Archive"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Base_Name",
                            "id": "#main/Metadata_Settings/Run_Base_Name"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Name",
                            "id": "#main/Metadata_Settings/Run_Name"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tag_Names",
                            "id": "#main/Metadata_Settings/Sample_Tag_Names"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/Metadata_Settings/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Start_Time/Start_Time",
                            "id": "#main/Metadata_Settings/Start_Time"
                        },
                        {
                            "source": "#main/Supplemental_Reference",
                            "id": "#main/Metadata_Settings/Supplemental_Reference"
                        },
                        {
                            "source": "#main/Targeted_Reference",
                            "id": "#main/Metadata_Settings/Targeted_Reference"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/Metadata_Settings/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#main/Metadata_Settings/Run_Metadata"
                    ],
                    "id": "#main/Metadata_Settings"
                },
                {
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "run": "#Metrics.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/GetDataTable/Cell_Type_Predictions",
                                "#main/ATAC_Compile_Results/Cell_Type_Predictions"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/Metrics/Cell_Type_Predictions"
                        },
                        {
                            "source": "#main/CalcGeneBodyCoverage/Gene_Body_Coverage_Json",
                            "id": "#main/Metrics/Gene_Body_Coverage_Json"
                        },
                        {
                            "source": [
                                "#main/ATAC_Compile_Results/Metrics_tar",
                                "#main/GetDataTable/Metrics_tar"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "first_non_null",
                            "id": "#main/Metrics/Metrics_tar"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Metrics/Run_Metadata"
                        },
                        {
                            "source": "#main/SampleTag/SampleTag_Metrics_Tar",
                            "id": "#main/Metrics/SampleTag_Metrics_Tar"
                        },
                        {
                            "source": "#main/CalcSeqSaturation/Seq_Saturation_Plot_Json",
                            "id": "#main/Metrics/Seq_Saturation_Plot_Json"
                        },
                        {
                            "source": "#main/VDJ_Compile_Results/vdjMetricsTar",
                            "id": "#main/Metrics/VDJ_Metrics_tar"
                        }
                    ],
                    "out": [
                        "#main/Metrics/Metrics_Summary",
                        "#main/Metrics/Metrics_Archive",
                        "#main/Metrics/Metrics_ATAC",
                        "#main/Metrics/Pipeline_Report_JSON",
                        "#main/Metrics/Pipeline_Report_HTML",
                        "#main/Metrics/output"
                    ],
                    "id": "#main/Metrics"
                },
                {
                    "label": "Miscellaneous Settings",
                    "run": "#MiscSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Custom_STAR_Params",
                            "id": "#main/Misc_Settings/_Custom_STAR_Params"
                        },
                        {
                            "source": "#main/Custom_bwa_mem2_Params",
                            "id": "#main/Misc_Settings/_Custom_bwa_mem2_Params"
                        },
                        {
                            "source": "#main/Exclude_Intronic_Reads",
                            "id": "#main/Misc_Settings/_Exclude_Intronic_Reads"
                        },
                        {
                            "source": "#main/Long_Reads",
                            "id": "#main/Misc_Settings/_Long_Reads"
                        },
                        {
                            "source": "#main/Multiple_Alignment_Gene_Selection",
                            "id": "#main/Misc_Settings/_Multiple_Alignment_Gene_Selection"
                        },
                        {
                            "source": "#main/Multiple_Alignment_Handling",
                            "id": "#main/Misc_Settings/_Multiple_Alignment_Handling"
                        }
                    ],
                    "out": [
                        "#main/Misc_Settings/Exclude_Intronic_Reads",
                        "#main/Misc_Settings/Long_Reads",
                        "#main/Misc_Settings/Custom_STAR_Params",
                        "#main/Misc_Settings/Custom_bwa_mem2_Params",
                        "#main/Misc_Settings/Multiple_Alignment_Handling",
                        "#main/Misc_Settings/Multiple_Alignment_Gene_Selection"
                    ],
                    "id": "#main/Misc_Settings"
                },
                {
                    "label": "Multiplexing Settings",
                    "run": "#MultiplexingSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Tag_Names",
                            "id": "#main/Multiplexing_Settings/_Sample_Tag_Names"
                        },
                        {
                            "source": "#main/Sample_Tags_Version",
                            "id": "#main/Multiplexing_Settings/_Sample_Tags_Version"
                        }
                    ],
                    "out": [
                        "#main/Multiplexing_Settings/Sample_Tag_Names",
                        "#main/Multiplexing_Settings/Sample_Tags_Version"
                    ],
                    "id": "#main/Multiplexing_Settings"
                },
                {
                    "label": "Name Settings",
                    "run": "#NameSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Reads_ATAC",
                            "valueFrom": "${ var fastqs = []; if(self) { for(var i = 0; i < self.length; i++) { fastqs.push(self[i].basename); } } return fastqs; }",
                            "id": "#main/Name_Settings/ATAC_Fastqs"
                        },
                        {
                            "source": "#main/Reads",
                            "valueFrom": "${ var fastqs = []; if(self) { for(var i = 0; i < self.length; i++) { fastqs.push(self[i].basename); } } return fastqs; }",
                            "id": "#main/Name_Settings/Bioproduct_Fastqs"
                        },
                        {
                            "source": "#main/Run_Name",
                            "id": "#main/Name_Settings/_Run_Name"
                        }
                    ],
                    "out": [
                        "#main/Name_Settings/Run_Name",
                        "#main/Name_Settings/Run_Base_Name"
                    ],
                    "id": "#main/Name_Settings"
                },
                {
                    "run": "#PeakAnnotation.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_ATAC",
                            "id": "#main/Peak_Annotation/Assay"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/Peak_Annotation/Gtf"
                        },
                        {
                            "source": "#main/QualCLAlign_ATAC/Peaks",
                            "id": "#main/Peak_Annotation/Peaks_bed"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Peak_Annotation/Run_Metadata"
                        }
                    ],
                    "out": [
                        "#main/Peak_Annotation/Peak_Annotation_TSV",
                        "#main/Peak_Annotation/output"
                    ],
                    "when": "$(inputs.Assay == \"ATAC\")",
                    "id": "#main/Peak_Annotation"
                },
                {
                    "label": "Putative Cell Calling Settings",
                    "run": "#PutativeCellSettings.cwl",
                    "in": [
                        {
                            "source": "#main/Cell_Calling_ATAC_Algorithm",
                            "id": "#main/Putative_Cell_Calling_Settings/_Cell_Calling_ATAC_Algorithm"
                        },
                        {
                            "source": "#main/Cell_Calling_Bioproduct_Algorithm",
                            "id": "#main/Putative_Cell_Calling_Settings/_Cell_Calling_Bioproduct_Algorithm"
                        },
                        {
                            "source": "#main/Cell_Calling_Data",
                            "id": "#main/Putative_Cell_Calling_Settings/_Cell_Calling_Data"
                        },
                        {
                            "source": "#main/Exact_Cell_Count",
                            "id": "#main/Putative_Cell_Calling_Settings/_Exact_Cell_Count"
                        },
                        {
                            "source": "#main/Expected_Cell_Count",
                            "id": "#main/Putative_Cell_Calling_Settings/_Expected_Cell_Count"
                        },
                        {
                            "source": "#main/Reads_ATAC",
                            "valueFrom": "${ if (self && self.length > 0) { return true; } else { return false; } }",
                            "id": "#main/Putative_Cell_Calling_Settings/_Reads_ATAC_Present"
                        },
                        {
                            "source": "#main/Reads",
                            "valueFrom": "${ if (self && self.length > 0) { return true; } else { return false; } }",
                            "id": "#main/Putative_Cell_Calling_Settings/_Reads_Present"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/Putative_Cell_Calling_Settings/_Sample_Tags_Version"
                        }
                    ],
                    "out": [
                        "#main/Putative_Cell_Calling_Settings/Cell_Calling_ATAC_Algorithm",
                        "#main/Putative_Cell_Calling_Settings/Cell_Calling_Bioproduct_Algorithm",
                        "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                        "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                        "#main/Putative_Cell_Calling_Settings/Expected_Cell_Count"
                    ],
                    "id": "#main/Putative_Cell_Calling_Settings"
                },
                {
                    "run": "#QualCLAlign.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_ATAC",
                            "id": "#main/QualCLAlign_ATAC/Assay"
                        },
                        {
                            "source": "#main/Bam_Settings/Generate_Bam",
                            "id": "#main/QualCLAlign_ATAC/Generate_Bam"
                        },
                        {
                            "source": "#main/CheckReference/Index",
                            "id": "#main/QualCLAlign_ATAC/Index"
                        },
                        {
                            "source": "#main/CheckReference/Checked_Predefined_Peaks",
                            "id": "#main/QualCLAlign_ATAC/Predefined_ATAC_Peaks"
                        },
                        {
                            "source": "#main/Reads_ATAC",
                            "id": "#main/QualCLAlign_ATAC/Reads"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Base_Name",
                            "id": "#main/QualCLAlign_ATAC/Run_Base_Name"
                        },
                        {
                            "source": [
                                "#main/Maximum_Threads",
                                "#main/GetMachineResources/Total_CPUs_Avail",
                                "#main/Assay_Settings/Assay_RNA"
                            ],
                            "valueFrom": "${ var [max_threads, total_cpus_avail, assay_rna] = self;\nif (max_threads != null) { return max_threads; }\nif (assay_rna == null) { return total_cpus_avail; }\nif (total_cpus_avail >= 48) { return parseInt(total_cpus_avail / 2); } else { return total_cpus_avail; } }",
                            "id": "#main/QualCLAlign_ATAC/Threads"
                        },
                        {
                            "source": "#main/Write_Filtered_Reads",
                            "id": "#main/QualCLAlign_ATAC/Write_Filtered_Reads"
                        },
                        {
                            "source": "#main/Misc_Settings/Custom_bwa_mem2_Params",
                            "id": "#main/QualCLAlign_ATAC/bwa_mem2_Params"
                        }
                    ],
                    "out": [
                        "#main/QualCLAlign_ATAC/Bead_Version",
                        "#main/QualCLAlign_ATAC/Libraries",
                        "#main/QualCLAlign_ATAC/ReadsList",
                        "#main/QualCLAlign_ATAC/BAMFiles",
                        "#main/QualCLAlign_ATAC/Fragments",
                        "#main/QualCLAlign_ATAC/Fragments_Index",
                        "#main/QualCLAlign_ATAC/Transposase_Sites",
                        "#main/QualCLAlign_ATAC/Transposase_Sites_Index",
                        "#main/QualCLAlign_ATAC/Peaks",
                        "#main/QualCLAlign_ATAC/Peaks_Index",
                        "#main/QualCLAlign_ATAC/QualCLAlignMetrics",
                        "#main/QualCLAlign_ATAC/UnifiedMetrics",
                        "#main/QualCLAlign_ATAC/Logs",
                        "#main/QualCLAlign_ATAC/Failed_Reads_CSVs",
                        "#main/QualCLAlign_ATAC/Reference_Genome_Size",
                        "#main/QualCLAlign_ATAC/Reference_Genome_Contig_Lengths"
                    ],
                    "when": "$(inputs.Assay == \"ATAC\")",
                    "id": "#main/QualCLAlign_ATAC"
                },
                {
                    "run": "#QualCLAlign.cwl",
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/QualCLAlign_RNA/Assay"
                        },
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/QualCLAlign_RNA/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/Index",
                            "id": "#main/QualCLAlign_RNA/Index"
                        },
                        {
                            "source": "#main/Misc_Settings/Long_Reads",
                            "id": "#main/QualCLAlign_RNA/Long_Reads"
                        },
                        {
                            "source": "#main/Reads",
                            "id": "#main/QualCLAlign_RNA/Reads"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Base_Name",
                            "id": "#main/QualCLAlign_RNA/Run_Base_Name"
                        },
                        {
                            "source": "#main/Misc_Settings/Custom_STAR_Params",
                            "id": "#main/QualCLAlign_RNA/STAR_Params"
                        },
                        {
                            "source": [
                                "#main/Maximum_Threads",
                                "#main/GetMachineResources/Total_CPUs_Avail",
                                "#main/Assay_Settings/Assay_ATAC"
                            ],
                            "valueFrom": "${ var [max_threads, total_cpus_avail, assay_atac] = self;\nif (max_threads != null) { return max_threads; }\nif (assay_atac == null) { return total_cpus_avail; }\nif (total_cpus_avail >= 48) { return parseInt(total_cpus_avail / 2); } else { return total_cpus_avail; } }",
                            "id": "#main/QualCLAlign_RNA/Threads"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/QualCLAlign_RNA/VDJ_Version"
                        },
                        {
                            "source": "#main/Write_Filtered_Reads",
                            "id": "#main/QualCLAlign_RNA/Write_Filtered_Reads"
                        }
                    ],
                    "out": [
                        "#main/QualCLAlign_RNA/Bead_Version",
                        "#main/QualCLAlign_RNA/Libraries",
                        "#main/QualCLAlign_RNA/ReadsList",
                        "#main/QualCLAlign_RNA/BAMFiles",
                        "#main/QualCLAlign_RNA/QualCLAlignMetrics",
                        "#main/QualCLAlign_RNA/Logs",
                        "#main/QualCLAlign_RNA/Failed_Reads_CSVs"
                    ],
                    "when": "$(inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\")",
                    "id": "#main/QualCLAlign_RNA"
                },
                {
                    "run": "#SampleTag.cwl",
                    "in": [
                        {
                            "source": "#main/AnnotateMolecules/Cell_Biop_Summary_List",
                            "pickValue": "all_non_null",
                            "id": "#main/SampleTag/Cell_Biop_Summary_List"
                        },
                        {
                            "source": "#main/CheckReference/Full_Genes",
                            "id": "#main/SampleTag/Full_Genes"
                        },
                        {
                            "source": [
                                "#main/AlignmentAnalysis/Gene_Name_ID_Json",
                                "#main/ATAC_Cell_by_Peak/Gene_Name_ID_Json"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/SampleTag/Gene_Name_ID_Json"
                        },
                        {
                            "source": [
                                "#main/GetDataTable/Cell_Order",
                                "#main/VDJ_Compile_Results/Cell_Order"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/SampleTag/Putative_Cells"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/SampleTag/Run_Metadata"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/SampleTag/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_bioproducts",
                            "id": "#main/SampleTag/num_bioproducts"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_cell_estimate",
                            "id": "#main/SampleTag/num_cell_estimate"
                        }
                    ],
                    "out": [
                        "#main/SampleTag/SampleTag_Metrics_Tar",
                        "#main/SampleTag/SampleTag_Calls",
                        "#main/SampleTag/SampleTag_CSVs",
                        "#main/SampleTag/SampleTag_perTagZips"
                    ],
                    "when": "$(inputs.Sample_Tags_Version != null)",
                    "id": "#main/SampleTag"
                },
                {
                    "run": {
                        "cwlVersion": "v1.2",
                        "class": "ExpressionTool",
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "inputs": [],
                        "outputs": [
                            {
                                "type": "string",
                                "id": "#main/Start_Time/run/Start_Time"
                            }
                        ],
                        "expression": "${\n  var today = new Date();\n  var date = today.toDateString();\n  var time = today.toLocaleTimeString('en-us', {timeZoneName: 'short'});\n  return ({Start_Time: date + ' ' + time});\n} "
                    },
                    "in": [],
                    "out": [
                        "#main/Start_Time/Start_Time"
                    ],
                    "id": "#main/Start_Time"
                },
                {
                    "run": {
                        "class": "ExpressionTool",
                        "cwlVersion": "v1.2",
                        "hints": [
                            {
                                "class": "ResourceRequirement",
                                "coresMin": 1,
                                "ramMin": 1024
                            }
                        ],
                        "inputs": [
                            {
                                "type": "string",
                                "id": "#main/SumMoleculesCorrected/run/Assay"
                            },
                            {
                                "type": {
                                    "type": "array",
                                    "items": "int"
                                },
                                "id": "#main/SumMoleculesCorrected/run/values"
                            }
                        ],
                        "outputs": [
                            {
                                "type": "int",
                                "id": "#main/SumMoleculesCorrected/run/Num_Mols_Corrected"
                            }
                        ],
                        "expression": "$({\n  Num_Mols_Corrected: inputs.values.reduce((a, b) => a + b, 0)\n})\n"
                    },
                    "in": [
                        {
                            "source": "#main/Assay_Settings/Assay_RNA",
                            "id": "#main/SumMoleculesCorrected/Assay"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Cell_Calling_Data",
                            "id": "#main/SumMoleculesCorrected/Cell_Calling_Data"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Num_Molecules_Corrected",
                            "id": "#main/SumMoleculesCorrected/values"
                        }
                    ],
                    "out": [
                        "#main/SumMoleculesCorrected/Num_Mols_Corrected"
                    ],
                    "when": "$( (inputs.Assay == \"WTA\" || inputs.Assay == \"Targeted\" ) && inputs.Cell_Calling_Data != 6 )",
                    "id": "#main/SumMoleculesCorrected"
                },
                {
                    "run": "#VDJ_Analyze_Reads_IG.cwl",
                    "when": "$(inputs.VDJ_Version != null && inputs.VDJ_Version != \"humanTCR\" && inputs.VDJ_Version != \"mouseTCR\")",
                    "in": [
                        {
                            "source": "#main/Maximum_Threads",
                            "id": "#main/VDJ_Analyze_Reads_IG/Maximum_Threads"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_valid_ig_reads",
                            "id": "#main/VDJ_Analyze_Reads_IG/Num_Valid_Reads_IG"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Analyze_Reads_IG/VDJ_Version"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/validIgReads",
                            "id": "#main/VDJ_Analyze_Reads_IG/Valid_Reads_Fastq_IG"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Analyze_Reads_IG/gatheredCalls"
                    ],
                    "id": "#main/VDJ_Analyze_Reads_IG"
                },
                {
                    "run": "#VDJ_Analyze_Reads_TCR.cwl",
                    "when": "$(inputs.VDJ_Version != null && inputs.VDJ_Version != \"humanBCR\" && inputs.VDJ_Version != \"mouseBCR\")",
                    "in": [
                        {
                            "source": "#main/Maximum_Threads",
                            "id": "#main/VDJ_Analyze_Reads_TCR/Maximum_Threads"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/num_valid_tcr_reads",
                            "id": "#main/VDJ_Analyze_Reads_TCR/Num_Valid_Reads_TCR"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Analyze_Reads_TCR/VDJ_Version"
                        },
                        {
                            "source": "#main/AlignmentAnalysis/validTcrReads",
                            "id": "#main/VDJ_Analyze_Reads_TCR/Valid_Reads_Fastq_TCR"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Analyze_Reads_TCR/gatheredCalls"
                    ],
                    "id": "#main/VDJ_Analyze_Reads_TCR"
                },
                {
                    "run": "#VDJ_Compile_Results.cwl",
                    "when": "$(inputs.VDJ_Version != null)",
                    "in": [
                        {
                            "source": "#main/AlignmentAnalysis/Seq_Metrics",
                            "id": "#main/VDJ_Compile_Results/Seq_Metrics"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Compile_Results/VDJ_Version"
                        },
                        {
                            "source": [
                                "#main/GetDataTable/Cell_Type_Predictions",
                                "#main/ATAC_Compile_Results/Cell_Type_Predictions"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/VDJ_Compile_Results/cellTypeMapping"
                        },
                        {
                            "valueFrom": "$([])",
                            "id": "#main/VDJ_Compile_Results/chainsToIgnore"
                        },
                        {
                            "source": "#main/Internal_Settings/VDJ_JGene_Evalue",
                            "id": "#main/VDJ_Compile_Results/evalueJgene"
                        },
                        {
                            "source": "#main/Internal_Settings/VDJ_VGene_Evalue",
                            "id": "#main/VDJ_Compile_Results/evalueVgene"
                        },
                        {
                            "source": "#main/VDJ_Analyze_Reads_IG/gatheredCalls",
                            "id": "#main/VDJ_Compile_Results/igCalls"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/VDJ_Compile_Results/metadata"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/VDJ_Compile_Results/putativeCells"
                        },
                        {
                            "source": "#main/VDJ_Analyze_Reads_TCR/gatheredCalls",
                            "id": "#main/VDJ_Compile_Results/tcrCalls"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Compile_Results/vdjCellsDatatable",
                        "#main/VDJ_Compile_Results/vdjCellsDatatableUncorrected",
                        "#main/VDJ_Compile_Results/vdjDominantContigsAIRR",
                        "#main/VDJ_Compile_Results/vdjUnfilteredContigsAIRR",
                        "#main/VDJ_Compile_Results/vdjAnnData",
                        "#main/VDJ_Compile_Results/vdjMetricsCsv",
                        "#main/VDJ_Compile_Results/vdjDbecFilterImages",
                        "#main/VDJ_Compile_Results/vdjMetricsTar",
                        "#main/VDJ_Compile_Results/Cell_Order"
                    ],
                    "id": "#main/VDJ_Compile_Results"
                },
                {
                    "label": "VDJ Settings",
                    "run": "#VDJ_Settings.cwl",
                    "in": [
                        {
                            "source": "#main/VDJ_Version",
                            "id": "#main/VDJ_Settings/_VDJ_Version"
                        }
                    ],
                    "out": [
                        "#main/VDJ_Settings/VDJ_Version"
                    ],
                    "id": "#main/VDJ_Settings"
                },
                {
                    "run": "#Version.cwl",
                    "in": [],
                    "out": [
                        "#main/Version/version"
                    ],
                    "id": "#main/Version"
                }
            ],
            "id": "#main"
        },
        {
            "requirements": [
                {
                    "listing": [
                        {
                            "entryname": "bam_files.txt",
                            "entry": "${\n    function getBamsInSortedOrder(inputs) {\n        // Create an associative array to hold the mapping from basename to full path\n        var fileMap = {};\n\n        // Extract basenames and map them to their full paths\n        for (var i = 0; i < inputs.BamFiles.length; i++) {\n            var file = inputs.BamFiles[i].path;\n            var basename = file.split('/').pop();\n            fileMap[basename] = file;\n        }\n\n        // Sort the basenames numerically\n        // This works because all bams share the same prefix and have a numerical id :\n        // i.e. foobar.0001.tagged.bam, foobar.0002.tagged.bam, foobar.0010.tagged.bam, etc.\n        //  so will be sorted by numerical ids\n        var sortedBasenames = Object.keys(fileMap).sort(function(a, b) {\n            return a.localeCompare(b, undefined, { numeric: true });\n        });\n\n        // Reconstruct the sorted full paths\n        var sortedBamFiles = sortedBasenames.map(function(basename) {\n            return fileMap[basename];\n        });\n\n        // Create a file of file names - 1 per line\n        return sortedBamFiles.join('\\n');\n    }\n\n    // For ATAC we cat the bams so need them in a particular order\n    if (inputs.Assay == \"ATAC\") {\n        return getBamsInSortedOrder(inputs);\n    }\n    else {\n        return inputs.BamFiles.map(function(file) {\n            return file.path;\n        }).join('\\n');\n    }\n}",
                            "writable": false
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "coresMin": 8,
                    "ramMin": 16000
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "MergeBam.sh"
            ],
            "stderr": "merge_bam.log",
            "inputs": [
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#MergeBAM.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#MergeBAM.cwl/BamFiles"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#MergeBAM.cwl/Generate_Bam"
                },
                {
                    "type": "File",
                    "loadContents": true,
                    "id": "#MergeBAM.cwl/Run_Metadata"
                }
            ],
            "arguments": [
                {
                    "position": 2,
                    "valueFrom": "${\n    var st_version = JSON.parse(inputs.Run_Metadata.contents).Sample_Tags_Version\n    if (st_version)\n    {\n      return st_version\n    } else\n    {\n      return \"None\"\n    }\n}"
                },
                {
                    "position": 3,
                    "valueFrom": "$(JSON.parse(inputs.Run_Metadata.contents).Run_Base_Name)"
                },
                {
                    "position": 4,
                    "valueFrom": "$(runtime.cores)"
                },
                {
                    "position": 5,
                    "valueFrom": "bam_files.txt"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "${ return \"*\" + JSON.parse(inputs.Run_Metadata.contents).Run_Base_Name + \"*.bam\" }"
                    },
                    "id": "#MergeBAM.cwl/Bam"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "${ return \"*\" + JSON.parse(inputs.Run_Metadata.contents).Run_Base_Name + \"*.bam.*\" }"
                    },
                    "id": "#MergeBAM.cwl/BamIndex"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#MergeBAM.cwl/log"
                }
            ],
            "id": "#MergeBAM.cwl"
        },
        {
            "class": "CommandLineTool",
            "baseCommand": "echo",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#Metadata.cwl/AbSeq_Reference"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "id": "#Metadata.cwl/Assay"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#Metadata.cwl/Bead_Version/Library",
                                    "type": "string"
                                },
                                {
                                    "name": "#Metadata.cwl/Bead_Version/bead_version",
                                    "type": "string"
                                }
                            ]
                        }
                    },
                    "id": "#Metadata.cwl/Bead_Version"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#Metadata.cwl/Cell_Calling_ATAC_Algorithm"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#Metadata.cwl/Cell_Calling_Bioproduct_Algorithm"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#Metadata.cwl/Cell_Calling_Data"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Custom_STAR_Params"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Custom_bwa_mem2_Params"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#Metadata.cwl/Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#Metadata.cwl/Exclude_Intronic_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#Metadata.cwl/Expected_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#Metadata.cwl/Generate_Bam"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#Metadata.cwl/Libraries"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#Metadata.cwl/Libraries_ATAC"
                },
                {
                    "type": "string",
                    "id": "#Metadata.cwl/Long_Reads"
                },
                {
                    "type": "string",
                    "id": "#Metadata.cwl/Pipeline_Name"
                },
                {
                    "type": "string",
                    "id": "#Metadata.cwl/Pipeline_Version"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Metadata.cwl/Predefined_ATAC_Peaks"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#Metadata.cwl/Reads"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#Metadata.cwl/Reads_ATAC"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Reference_Archive"
                },
                {
                    "type": "string",
                    "id": "#Metadata.cwl/Run_Base_Name"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Run_Name"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#Metadata.cwl/Sample_Tag_Names"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Start_Time"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#Metadata.cwl/Supplemental_Reference"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#Metadata.cwl/Targeted_Reference"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "run_metadata.json"
                    },
                    "id": "#Metadata.cwl/Run_Metadata"
                }
            ],
            "stdout": "run_metadata.json",
            "arguments": [
                {
                    "prefix": ""
                },
                {
                    "shellQuote": true,
                    "valueFrom": "${\n  var metadata = inputs;\n  var all_bv = {};\n  var customer_bv = \"Original (V1)\";\n  var detected_bv = \"V1\";\n  for (var i = 0; i < inputs.Bead_Version.length; i++) {\n      var BeadVer = inputs.Bead_Version[i];\n      var Library = BeadVer[\"Library\"];\n      var bead_version = BeadVer[\"bead_version\"];\n      all_bv[Library] = bead_version;\n      var short_bv =  bead_version.substring(0, 5);\n      if (short_bv == \"Enh\") {\n        customer_bv = \"Enhanced\";\n        detected_bv = \"Enh\";\n      }\n      else if (short_bv == \"EnhV2\") {\n        customer_bv = \"Enhanced V2/V3\";\n        detected_bv = \"EnhV2\";\n      }\n  }\n  metadata[\"Bead_Version\"] = all_bv;\n  metadata[\"Bead_Version_Detected\"] = detected_bv;\n\n  var pipeline_name = inputs.Pipeline_Name;\n  var version = inputs.Pipeline_Version;\n  var time = inputs.Start_Time;\n  var libraries = inputs.Libraries;\n  if(libraries == null){\n      libraries = [\"None\"];\n  }\n  var libraries_atac = inputs.Libraries_ATAC;\n  if(libraries_atac == null){\n      libraries_atac = [\"None\"];\n  }\n\n  var i = 0;\n  var refs_mrna = [];\n  var mrna_name = \"mRNA Reference\";\n  if (inputs.Targeted_Reference != null) {\n      mrna_name = \"Targeted Reference\";\n      for (i = 0; i < inputs.Targeted_Reference.length; i++) {\n          if (inputs.Targeted_Reference[i] != null) {\n              refs_mrna.push(inputs.Targeted_Reference[i][\"basename\"]);\n          }\n      }\n  }\n  if(inputs.Reference_Archive != null){\n      mrna_name = \"Reference Archive\";\n      refs_mrna.push(inputs.Reference_Archive);\n  }\n  if (refs_mrna.length == 0) {\n      refs_mrna = [\"None\"];\n  }\n\n  var refs_abseq_inputs = [];\n  if (inputs.AbSeq_Reference != null) {\n      refs_abseq_inputs = refs_abseq_inputs.concat(inputs.AbSeq_Reference);\n  }\n  var refs_abseq = [];\n  if (refs_abseq_inputs.length > 0) {\n      for (i = 0; i < refs_abseq_inputs.length; i++) {\n          if (refs_abseq_inputs[i] != null) {\n              refs_abseq.push(refs_abseq_inputs[i][\"basename\"]);\n          }\n      }\n  }\n  else {\n      refs_abseq = [\"None\"];\n  }\n\n  var refs_supp_inputs = [];\n  if (inputs.Supplemental_Reference != null) {\n      refs_supp_inputs = refs_supp_inputs.concat(inputs.AbSeq_Reference);\n  }\n  var refs_supp = [];\n  if (refs_supp_inputs.length > 0) {\n      for (i = 0; i < refs_supp_inputs.length; i++) {\n          if (refs_supp_inputs[i] != null) {\n              refs_supp.push(refs_supp_inputs[i][\"basename\"]);\n          }\n      }\n  }\n  else {\n      refs_supp = [\"None\"];\n  }\n\n  if (inputs.Predefined_ATAC_Peaks != null) {\n      var predef_atac_peaks = inputs.Predefined_ATAC_Peaks[\"basename\"];\n  } else {\n      var predef_atac_peaks = \"None\";\n  }\n\n  var parameters = [];\n  if(inputs.Sample_Tags_Version != null){\n      var tags = \"Sample Tag Version: \" + inputs.Sample_Tags_Version;\n  } else{ \n      var tags = \"Sample Tag Version: None\";\n  }\n  parameters.push(tags);\n\n  if(inputs.Sample_Tag_Names != null){\n      var tag_names = inputs.Sample_Tag_Names.join(\"; \");\n      var tag_list = \"Sample Tag Names: \" + tag_names;\n  } else{\n      var tag_list = \"Sample Tag Names: None\";\n  }\n  parameters.push(tag_list);\n \n  if(inputs.VDJ_Version != null){\n      var vdj = \"VDJ Version: \" + inputs.VDJ_Version;\n  } else{ \n      var vdj = \"VDJ Version: None\";\n  }\n  parameters.push(vdj);\n\n  if (inputs.Cell_Calling_Data == 0) {\n      var call = \"Putative Cell Calling Data: mRNA\";\n  } else if (inputs.Cell_Calling_Data == 1) {\n      var call = \"Putative Cell Calling Data: AbSeq\";\n  } else if (inputs.Cell_Calling_Data == 2) {\n      var call = \"Putative Cell Calling Data: mRNA_and_AbSeq\";\n  } else if (inputs.Cell_Calling_Data == 3) {\n      var call = \"Putative Cell Calling Data: mRNA_and_ATAC\";\n  } else if (inputs.Cell_Calling_Data == 4) {\n      var call = \"Putative Cell Calling Data: AbSeq_and_ATAC\";\n  } else if (inputs.Cell_Calling_Data == 5) {\n      var call = \"Putative Cell Calling Data: ATAC\";\n  } else if (inputs.Cell_Calling_Data == 6) {\n      var call = \"Putative_Cell_Calling_Data: VDJ\";\n  } else {\n      var call = \"Putative Cell Calling Data: None\";\n  }\n  parameters.push(call);\n\n  if (inputs.Cell_Calling_Bioproduct_Algorithm != null) {\n      var bioproduct_alg = \"Bioproduct Cell Calling Algorithm: \" + inputs.Cell_Calling_Bioproduct_Algorithm;\n  } else {\n      var bioproduct_alg = \"Bioproduct Cell Calling Algorithm: None\";\n  }\n  parameters.push(bioproduct_alg);\n\n  if (inputs.Cell_Calling_ATAC_Algorithm != null) {\n      var atac_alg = \"ATAC Cell Calling Algorithm: \" + inputs.Cell_Calling_ATAC_Algorithm;\n  } else {\n      var atac_alg = \"ATAC Cell Calling Algorithm: None\";\n  }\n  parameters.push(atac_alg);\n\n  if(inputs.Exclude_Intronic_Reads){\n      var introns = \"Exclude Intronic Reads: On\";\n  } else{\n      var introns = \"Exclude Intronic Reads: Off\";\n  }\n  parameters.push(introns);\n\n  if(inputs.Generate_Bam){\n      var generateBam = \"Generate Bam: On\";\n  } else{\n      var generateBam = \"Generate Bam: Off\";\n  }\n  parameters.push(generateBam);\n\n  if(inputs.Exact_Cell_Count != null){\n      var exactCells = \"Exact Cell Count: \" + inputs.Exact_Cell_Count;\n  } else{\n      var exactCells = \"Exact Cell Count: None\";\n  }\n  parameters.push(exactCells);\n\n  if(inputs.Expected_Cell_Count != null){\n      var expectedCells = \"Expected Cell Count: \" + inputs.Expected_Cell_Count;\n  } else{\n      var expectedCells = \"Expected Cell Count: None\";\n  }\n  parameters.push(expectedCells);\n\n  var longReads = \"Long Reads: \" + inputs.Long_Reads;\n  parameters.push(longReads);\n\n  if (inputs.Custom_STAR_Params != null)\n  {\n    var starParams = \"Custom STAR Params: \" + inputs.Custom_STAR_Params;\n  } else {\n    var starParams = \"Custom STAR Params: None\"; \n  }\n  parameters.push(starParams);\n\n  if (inputs.Custom_bwa_mem2_Params != null)\n  {\n    var bwaParams = \"Custom bwa-mem2 Params: \" + inputs.Custom_bwa_mem2_Params;\n  } else {\n    var bwaParams = \"Custom bwa-mem2 Params: None\";\n  }\n  parameters.push(bwaParams);\n\n  var run_name = inputs.Run_Name;\n  var run_base_name = inputs.Run_Base_Name;\n\n  var header = [\"####################\"];\n  header.push(\"## \" + pipeline_name + \" Version \" + version);\n  header.push(\"## Analysis Date - \" + time);\n  header.push(\"## Libraries - Bioproduct Libraries: \" + libraries.join('; ') + \" | ATAC Libraries: \" + libraries_atac.join('; ') + \" | Bead version detected: \" + customer_bv);\n  header.push(\"## References - \" + mrna_name + \": \" + refs_mrna.join('; ') + \" | AbSeq Reference: \" + refs_abseq.join('; ') + \" | Supplemental Reference: \" + refs_supp.join('; ') + \" | ATAC Predefined Peak Regions: \" + predef_atac_peaks);\n  header.push(\"## Parameters - \" + parameters.join(' | '));\n  header.push(\"####################\");\n  metadata[\"Output_Header\"] = header;\n  metadata[\"Run_Name\"] = run_name;   \n  metadata[\"Run_Base_Name\"] = run_base_name;\n\n  var metadata_json = JSON.stringify(metadata, null, 2);\n\n  return metadata_json;\n}\n"
                }
            ],
            "id": "#Metadata.cwl"
        },
        {
            "hints": [
                {
                    "dockerPull": "bdgenomics/rhapsody:4.0",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_metrics.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "prefix": "--cell-type-file",
                        "itemSeparator": ","
                    },
                    "id": "#Metrics.cwl/Cell_Type_Predictions"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--gene-body-coverage-json"
                    },
                    "id": "#Metrics.cwl/Gene_Body_Coverage_Json"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--metrics-tar"
                    },
                    "id": "#Metrics.cwl/Metrics_tar"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#Metrics.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--sample-tag-metrics-tar"
                    },
                    "id": "#Metrics.cwl/SampleTag_Metrics_Tar"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--seq-saturation-plot-json"
                    },
                    "id": "#Metrics.cwl/Seq_Saturation_Plot_Json"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-metrics-tar"
                    },
                    "id": "#Metrics.cwl/VDJ_Metrics_tar"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": [
                            "*_ATAC_Metrics.json",
                            "*_ATAC_Metrics.csv"
                        ]
                    },
                    "id": "#Metrics.cwl/Metrics_ATAC"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "internal-metrics-archive.tar.gz"
                    },
                    "id": "#Metrics.cwl/Metrics_Archive"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Metrics_Summary.csv"
                    },
                    "id": "#Metrics.cwl/Metrics_Summary"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Pipeline_Report.html"
                    },
                    "id": "#Metrics.cwl/Pipeline_Report_HTML"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_Pipeline_Report.json"
                    },
                    "id": "#Metrics.cwl/Pipeline_Report_JSON"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#Metrics.cwl/output"
                }
            ],
            "id": "#Metrics.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MiscSettings.cwl/_Custom_STAR_Params"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MiscSettings.cwl/_Custom_bwa_mem2_Params"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#MiscSettings.cwl/_Exclude_Intronic_Reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#MiscSettings.cwl/_Long_Reads"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#MiscSettings.cwl/_Multiple_Alignment_Gene_Selection"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#MiscSettings.cwl/_Multiple_Alignment_Handling"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MiscSettings.cwl/Custom_STAR_Params"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MiscSettings.cwl/Custom_bwa_mem2_Params"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#MiscSettings.cwl/Exclude_Intronic_Reads"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MiscSettings.cwl/Long_Reads"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#MiscSettings.cwl/Multiple_Alignment_Gene_Selection"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#MiscSettings.cwl/Multiple_Alignment_Handling"
                }
            ],
            "expression": "${\n  // the exclude intronic reads flag defaults to false\n  var excludeIntronicReads = false;\n  // the user can set the flag to exclude intronic reads\n  if (inputs._Exclude_Intronic_Reads) {\n    excludeIntronicReads = inputs._Exclude_Intronic_Reads;\n  }\n\n  // Use_Long_Reads default is autodetect, which happens in CheckFastqs\n  // User can set this explicitly true or false\n  // Convert boolean results to string.  null -> \"auto\", true -> \"true\", false -> \"false\"\n  var longReads = \"auto\";\n  if (inputs._Long_Reads !== null) {\n    if (inputs._Long_Reads) {\n      longReads = \"true\";\n    }\n    else {\n      longReads = \"false\";\n    }\n  }\n\n  return ({\n    Exclude_Intronic_Reads: excludeIntronicReads,\n    Long_Reads: longReads,\n    Custom_STAR_Params: inputs._Custom_STAR_Params,\n    Custom_bwa_mem2_Params: inputs._Custom_bwa_mem2_Params,\n    Multiple_Alignment_Handling: inputs._Multiple_Alignment_Handling,\n    Multiple_Alignment_Gene_Selection: inputs._Multiple_Alignment_Gene_Selection\n  });\n\n\n}",
            "id": "#MiscSettings.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "default": "Targeted",
                    "id": "#MultiplexingSettings.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/_Sample_Tag_Names"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#MultiplexingSettings.cwl/_Sample_Tags_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/Sample_Tag_Names"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MultiplexingSettings.cwl/Sample_Tags_Version"
                }
            ],
            "expression": "${\n  var enumifiedSampleTagsVersion = null;\n  if (inputs._Sample_Tags_Version) {\n  var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();\n\n  if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')\n  {\n    enumifiedSampleTagsVersion = 'hs';\n  }\n  else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')\n  {\n    enumifiedSampleTagsVersion = 'mm';\n  }\n  else if (_Sample_Tags_Version.indexOf('flex') >= 0)\n  {\n    enumifiedSampleTagsVersion = 'flex';\n  }\n  else if (_Sample_Tags_Version.indexOf('nuclei') >= 0)\n  {\n    if (_Sample_Tags_Version.indexOf('atac') >= 0)\n    {\n      enumifiedSampleTagsVersion = 'nuclei_atac_only';\n    }\n    else\n    {\n      enumifiedSampleTagsVersion = 'nuclei_includes_mrna';\n    }\n  }\n  else if (_Sample_Tags_Version === 'no multiplexing')\n  {\n    enumifiedSampleTagsVersion = null;\n  }\n  else\n  {\n    throw new Error(\"Cannot parse Sample Tag Version: \" + inputs._Sample_Tags_Version);\n  }\n  }\n  var newSampleTags = null;\n  var warning = '';\n  var malformedSampleTags = [];\n  if (inputs._Sample_Tag_Names) {\n    var listSampleTags = inputs._Sample_Tag_Names;\n    var newSampleTags = [];\n\n    // loop over the custom sample tags\n    for (var i = 0; i < listSampleTags.length; i++) {\n      var sampleTag = listSampleTags[i];\n      var newSampleTagItems = [];\n\n      // the custom sample tag must include a dash between the number and the name\n      if (sampleTag.includes(\"-\")) {\n          var sampleTagItemsList = sampleTag.split(\"-\");\n\n          // loop over the items separated by a dash\n          for (var j = 0; j < sampleTagItemsList.length; j++) {\n              var sampleTagItem = sampleTagItemsList[j]\n\n              // remove the leading and trailing spaces\n              newSampleTagItems.push(sampleTagItem.trim());\n          }\n          // after removing all leading and trailing spaces from the items separated by a dash\n          var newSampleTag = newSampleTagItems.join(\"-\");\n\n          // only alpha-numeric and + and - are allowed, replace everything else with underscore\n          var finalSampleTag = newSampleTag.replace(/[^A-Za-z0-9-+]/g,\"_\");\n          newSampleTags.push(finalSampleTag);\n      } else {\n          malformedSampleTags.push(sampleTag);\n      }\n    }\n  }\n\n  // add a warning about malformed sample tags\n  if (malformedSampleTags.length > 0) {\n      warning = \"The following sample tag names do not contain the required dash symbol (-) and will be skipped: [\" + malformedSampleTags.join(\",\") + \"]\";\n  }\n\n  return ({\n  Sample_Tag_Names: newSampleTags,\n  Sample_Tags_Version: enumifiedSampleTagsVersion,\n  Warning: warning\n  });\n}",
            "id": "#MultiplexingSettings.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "doc": "NameSettings sets the Run_Name variable that is used as a common prefix to name output files. If the user has specified a Run_Name, it is cleaned up or one is set based on the Bioproduct/ATAC fastq filenames.\n",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                }
            ],
            "baseCommand": [
                "python",
                "-c"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "default": [],
                    "id": "#NameSettings.cwl/ATAC_Fastqs"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "default": [],
                    "id": "#NameSettings.cwl/Bioproduct_Fastqs"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "",
                    "id": "#NameSettings.cwl/_Run_Name"
                }
            ],
            "arguments": [
                {
                    "position": 2,
                    "valueFrom": "import sys\nfrom mist.apps import CheckFastqs\nCheckFastqs.write_run_name_json(sys.argv[1], sys.argv[2], sys.argv[3])\n",
                    "shellQuote": true
                },
                {
                    "position": 3,
                    "valueFrom": "$(inputs._Run_Name)",
                    "shellQuote": true
                },
                {
                    "position": 4,
                    "valueFrom": "$(inputs.Bioproduct_Fastqs.join(\",\"))",
                    "shellQuote": true
                },
                {
                    "position": 5,
                    "valueFrom": "$(inputs.ATAC_Fastqs.join(\",\"))",
                    "shellQuote": true
                }
            ],
            "outputs": [
                {
                    "type": "string",
                    "outputBinding": {
                        "glob": "run_base_name.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents)['Run_Base_Name'])"
                    },
                    "id": "#NameSettings.cwl/Run_Base_Name"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "outputBinding": {
                        "outputEval": "$(inputs._Run_Name)"
                    },
                    "id": "#NameSettings.cwl/Run_Name"
                }
            ],
            "id": "#NameSettings.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 16000
                }
            ],
            "baseCommand": [
                "mist_peak_annotation.py"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "id": "#PeakAnnotation.cwl/Gtf"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--peaks_bed_file"
                    },
                    "secondaryFiles": {
                        "pattern": ".tbi",
                        "required": null
                    },
                    "id": "#PeakAnnotation.cwl/Peaks_bed"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#PeakAnnotation.cwl/Run_Metadata"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.tsv.gz"
                    },
                    "id": "#PeakAnnotation.cwl/Peak_Annotation_TSV"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "id": "#PeakAnnotation.cwl/output"
                }
            ],
            "id": "#PeakAnnotation.cwl"
        },
        {
            "class": "CommandLineTool",
            "baseCommand": "echo",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Cell_Calling_ATAC_Algorithm"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Cell_Calling_Bioproduct_Algorithm"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Cell_Calling_Data"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Expected_Cell_Count"
                },
                {
                    "type": "boolean",
                    "id": "#PutativeCellSettings.cwl/_Reads_ATAC_Present"
                },
                {
                    "type": "boolean",
                    "id": "#PutativeCellSettings.cwl/_Reads_Present"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Sample_Tags_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "outputBinding": {
                        "glob": "putative_cell_settings.json",
                        "loadContents": true,
                        "outputEval": "${ return JSON.parse(self[0].contents).Cell_Calling_ATAC_Algorithm; }"
                    },
                    "id": "#PutativeCellSettings.cwl/Cell_Calling_ATAC_Algorithm"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "outputBinding": {
                        "glob": "putative_cell_settings.json",
                        "loadContents": true,
                        "outputEval": "${ return JSON.parse(self[0].contents).Cell_Calling_Bioproduct_Algorithm; }"
                    },
                    "id": "#PutativeCellSettings.cwl/Cell_Calling_Bioproduct_Algorithm"
                },
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "putative_cell_settings.json",
                        "loadContents": true,
                        "outputEval": "${ return JSON.parse(self[0].contents).Cell_Calling_Data; }"
                    },
                    "id": "#PutativeCellSettings.cwl/Cell_Calling_Data"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "outputBinding": {
                        "glob": "putative_cell_settings.json",
                        "loadContents": true,
                        "outputEval": "${ return JSON.parse(self[0].contents).Exact_Cell_Count; }"
                    },
                    "id": "#PutativeCellSettings.cwl/Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "outputBinding": {
                        "glob": "putative_cell_settings.json",
                        "loadContents": true,
                        "outputEval": "${ return JSON.parse(self[0].contents).Expected_Cell_Count; }"
                    },
                    "id": "#PutativeCellSettings.cwl/Expected_Cell_Count"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "putative_cell_settings.json"
                    },
                    "id": "#PutativeCellSettings.cwl/PutativeCellSettings"
                }
            ],
            "stdout": "putative_cell_settings.json",
            "arguments": [
                {
                    "prefix": ""
                },
                {
                    "shellQuote": true,
                    "valueFrom": "${\n  var settings = inputs;\n  var errorMessage = \"No error\";\n\n  var cellCallingATACAlg = null;\n  // the default cell calling algorithm for ATAC is basic\n  if (inputs._Reads_ATAC_Present) {\n      cellCallingATACAlg = \"Basic\";\n  }\n  // the user can choose the ATAC cell calling algorithm\n  if (inputs._Cell_Calling_ATAC_Algorithm) {\n      cellCallingATACAlg = inputs._Cell_Calling_ATAC_Algorithm;\n  }\n\n  var cellCallingBioproductAlg = null;\n  // the default cell calling algorithm for bioproducts is basic\n  if (inputs._Reads_Present) {\n      cellCallingBioproductAlg = \"Basic\";\n  }\n  // the user can choose the bioproducts cell calling algorithm\n  if (inputs._Cell_Calling_Bioproduct_Algorithm) {\n      cellCallingBioproductAlg = inputs._Cell_Calling_Bioproduct_Algorithm;\n  }\n\n  // the default cell calling data depends on the data that is provided\n  // the overall default is mRNA data\n  var cellCallingDataInt = 0;\n  // if reads and ATAC reads are present\n  if (inputs._Reads_Present && inputs._Reads_ATAC_Present) {\n      // if the sample tag version is set to ATAC-only,\n      // then the reads are SMK reads, default to ATAC cell calling\n      if (inputs._Sample_Tags_Version && inputs._Sample_Tags_Version === \"nuclei_atac_only\") {\n          cellCallingDataInt = 5;\n      }\n      // otherwise default to joint cell calling\n      else {\n          cellCallingDataInt = 3;\n      }\n  }\n  // if no mRNA or AbSeq data is present, but ATAC data is, then default to ATAC\n  else if (!inputs._Reads_Present && inputs._Reads_ATAC_Present) {\n      cellCallingDataInt = 5;\n  }\n\n  // convert the Cell_Calling_Data from a string to an integer\n  if (inputs._Cell_Calling_Data) {\n    if (inputs._Cell_Calling_Data === \"mRNA\") {\n      if (!inputs._Reads_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"mRNA\" was selected but no mRNA Reads were provided.';\n      } else {\n        cellCallingDataInt = 0;\n      }\n    }\n    else if (inputs._Cell_Calling_Data === \"AbSeq\") {\n      if (!inputs._Reads_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"AbSeq\" was selected but no AbSeq Reads were provided.';\n      }\n      cellCallingDataInt = 1;\n    }\n    else if (inputs._Cell_Calling_Data === \"mRNA_and_AbSeq\") {\n      if (!inputs._Reads_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"mRNA_and_AbSeq\" was selected but no mRNA/AbSeq Reads were provided.';\n      }\n      cellCallingDataInt = 2;\n    }\n    else if (inputs._Cell_Calling_Data === \"mRNA_and_ATAC\") {\n      if (!inputs._Reads_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"mRNA_and_ATAC\" was selected but no mRNA Reads were provided.';\n      } else if (!inputs._Reads_ATAC_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"mRNA_and_ATAC\" was selected but no ATAC Reads were provided.';\n      }\n      cellCallingDataInt = 3;\n    }\n    else if (inputs._Cell_Calling_Data === \"AbSeq_and_ATAC\") {\n      if (!inputs._Reads_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"AbSeq_and_ATAC\" was selected but no AbSeq Reads were provided.';\n      } else if (!inputs._Reads_ATAC_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"AbSeq_and_ATAC\" was selected but no ATAC Reads were provided.';\n      }\n      cellCallingDataInt = 4;\n    }\n    else if (inputs._Cell_Calling_Data === \"ATAC\") {\n      if (!inputs._Reads_ATAC_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"ATAC\" was selected but no ATAC Reads were provided.';\n      }\n      cellCallingDataInt = 5;\n    }\n    else if (inputs._Cell_Calling_Data === \"VDJ\") {\n      if (!inputs._Reads_Present) {\n        errorMessage = 'The \"Cell Calling Data\" option \"VDJ\" was selected but no mRNA Reads were provided.';\n      }\n      cellCallingDataInt = 6;\n    }\n    // if reads and ATAC reads are present\n    if (inputs._Reads_Present && inputs._Reads_ATAC_Present) {\n      // if the sample tag version is set to ATAC-only,\n      // then the reads are probably SMK reads\n      // create error if ATAC is not selected for cell calling\n      if (inputs._Sample_Tags_Version && inputs._Sample_Tags_Version === \"nuclei_atac_only\" && inputs._Cell_Calling_Data != \"ATAC\") {\n        errorMessage = 'Incompatible input parameters detected: the \"Cell Calling Data\" option \"' + inputs._Cell_Calling_Data + '\" and the \"Sample Tag Version\" option \"Nuclei - ATAC-Only\" are incompatible.';\n      }\n    }\n  }\n  // check the exact cell count\n  if (inputs._Exact_Cell_Count) {\n    if (inputs._Exact_Cell_Count < 1) {\n      errorMessage = \"Exact cell count must be an integer greater than 0, value received: \" + inputs._Exact_Cell_Count;\n    }\n  }\n  // check if there is an error\n  if (errorMessage != \"No error\") {\n    // If there is an error, force CWL to show it:\n    //   - \"Cell_Calling_Data\" is a required output\n    //   - setting it to null will cause a CWL error\n    //   - the error message will be shown in the json\n    cellCallingDataInt = null;\n  }\n\n  settings[\"Cell_Calling_ATAC_Algorithm\"] = cellCallingATACAlg;\n  settings[\"Cell_Calling_Bioproduct_Algorithm\"] = cellCallingBioproductAlg;\n  settings[\"Cell_Calling_Data\"] = cellCallingDataInt;\n  settings[\"Expected_Cell_Count\"] = inputs._Expected_Cell_Count;\n  settings[\"Exact_Cell_Count\"] = inputs._Exact_Cell_Count;\n  settings[\"Error\"] = errorMessage;\n\n  var settings_json = JSON.stringify(settings, null, 2);\n\n  return settings_json;\n  }\n"
                }
            ],
            "id": "#PutativeCellSettings.cwl"
        },
        {
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "doc": "CheckFastqs does several quality control routines including: ensuring that read pair file names are formatted correctly and contain a read pair mate;  QualCLAlign stage of the Rain pipeline overlaps read pairs and then performs a series of filters and mappings to reduce valid reads into a single FastQ file to be fed into the aligner. The R2 reads are annotated with cell index and UMI information derived from the R1 read.\n",
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "coresMin": "${ if(inputs.Threads){ return inputs.Threads; } else{ return 8; } }",
                    "ramMin": 48000
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_run_qualclalign.py"
            ],
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--alignment-compression-threads"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#QualCLAlign.cwl/Alignment_Compression_threads"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--assay"
                    },
                    "id": "#QualCLAlign.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--bgzf-threads"
                    },
                    "id": "#QualCLAlign.cwl/BGZF_Threads"
                },
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualCLAlign.cwl/Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "prefix": "--split-atac-bam"
                    },
                    "id": "#QualCLAlign.cwl/Generate_Bam"
                },
                {
                    "inputBinding": {
                        "prefix": "--index"
                    },
                    "type": "Directory",
                    "id": "#QualCLAlign.cwl/Index"
                },
                {
                    "inputBinding": {
                        "prefix": "--use-star-long"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#QualCLAlign.cwl/Long_Reads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--predefined-peaks"
                    },
                    "secondaryFiles": {
                        "pattern": ".tbi",
                        "required": null
                    },
                    "id": "#QualCLAlign.cwl/Predefined_ATAC_Peaks"
                },
                {
                    "inputBinding": {
                        "prefix": "--reader-annotation-threads"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#QualCLAlign.cwl/Reader_Annotation_Threads"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "prefix": "--reads",
                        "itemSeparator": ","
                    },
                    "id": "#QualCLAlign.cwl/Reads"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--run-name"
                    },
                    "id": "#QualCLAlign.cwl/Run_Base_Name"
                },
                {
                    "inputBinding": {
                        "prefix": "--star-params",
                        "shellQuote": true
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#QualCLAlign.cwl/STAR_Params"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#QualCLAlign.cwl/Threads"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-version"
                    },
                    "id": "#QualCLAlign.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "prefix": "--write-filtered-read-pairs"
                    },
                    "id": "#QualCLAlign.cwl/Write_Filtered_Reads"
                },
                {
                    "inputBinding": {
                        "prefix": "--bwa-mem2-params",
                        "shellQuote": true
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#QualCLAlign.cwl/bwa_mem2_Params"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*.bam"
                    },
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#QualCLAlign.cwl/BAMFiles"
                },
                {
                    "type": {
                        "type": "array",
                        "items": {
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#QualCLAlign.cwl/Bead_Version/Library",
                                    "type": "string"
                                },
                                {
                                    "name": "#QualCLAlign.cwl/Bead_Version/bead_version",
                                    "type": "string"
                                }
                            ]
                        }
                    },
                    "outputBinding": {
                        "glob": "bead_version.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).BeadVersion)\n"
                    },
                    "id": "#QualCLAlign.cwl/Bead_Version"
                },
                {
                    "outputBinding": {
                        "glob": "*.failedReads.csv.gz"
                    },
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#QualCLAlign.cwl/Failed_Reads_CSVs"
                },
                {
                    "outputBinding": {
                        "glob": "fastq_read_pairs.json"
                    },
                    "type": "File",
                    "id": "#QualCLAlign.cwl/Fastq_read_pairs"
                },
                {
                    "outputBinding": {
                        "glob": "*_ATAC_Fragments.bed.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "secondaryFiles": [
                        {
                            "pattern": ".tbi",
                            "required": null
                        }
                    ],
                    "id": "#QualCLAlign.cwl/Fragments"
                },
                {
                    "outputBinding": {
                        "glob": "*_ATAC_Fragments.bed.gz.tbi"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualCLAlign.cwl/Fragments_Index"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "outputBinding": {
                        "glob": "bead_version.json",
                        "loadContents": true,
                        "outputEval": "${\n  var obj = JSON.parse(self[0].contents);\n  var libraries = [];\n  var beadLibs = obj.BeadVersion\n  for (var i in beadLibs){\n    if (libraries.indexOf(beadLibs[i][\"Library\"]) == -1){ \n      libraries.push(beadLibs[i][\"Library\"]);\n    }\n  }\n  libraries.sort();\n  return libraries\n}\n"
                    },
                    "id": "#QualCLAlign.cwl/Libraries"
                },
                {
                    "outputBinding": {
                        "glob": "*logs.tar.gz"
                    },
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#QualCLAlign.cwl/Logs"
                },
                {
                    "outputBinding": {
                        "glob": "*_ATAC_Peaks.bed.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "secondaryFiles": [
                        {
                            "pattern": ".tbi",
                            "required": null
                        }
                    ],
                    "id": "#QualCLAlign.cwl/Peaks"
                },
                {
                    "outputBinding": {
                        "glob": "*_ATAC_Peaks.bed.gz.tbi"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualCLAlign.cwl/Peaks_Index"
                },
                {
                    "outputBinding": {
                        "glob": "*ReadQualityMetrics.json"
                    },
                    "type": "File",
                    "id": "#QualCLAlign.cwl/QualCLAlignMetrics"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "outputBinding": {
                        "outputEval": "${  \n  var reads = []; \n  var files = inputs.Reads\n  for (var i in files){\n      reads.push(files[i][\"basename\"]);\n  }\n  reads.sort();\n  return(reads)\n}\n"
                    },
                    "id": "#QualCLAlign.cwl/ReadsList"
                },
                {
                    "outputBinding": {
                        "glob": "genome_contig_lengths.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualCLAlign.cwl/Reference_Genome_Contig_Lengths"
                },
                {
                    "type": "string",
                    "outputBinding": {
                        "glob": "genome_size.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return \"0\"; } return JSON.parse(self[0].contents); }"
                    },
                    "id": "#QualCLAlign.cwl/Reference_Genome_Size"
                },
                {
                    "outputBinding": {
                        "glob": "*_ATAC_Transposase_Sites.bed.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "secondaryFiles": [
                        {
                            "pattern": ".tbi",
                            "required": null
                        }
                    ],
                    "id": "#QualCLAlign.cwl/Transposase_Sites"
                },
                {
                    "outputBinding": {
                        "glob": "*_ATAC_Transposase_Sites.bed.gz.tbi"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualCLAlign.cwl/Transposase_Sites_Index"
                },
                {
                    "outputBinding": {
                        "glob": "*_ATAC_UnifiedMetrics.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualCLAlign.cwl/UnifiedMetrics"
                }
            ],
            "id": "#QualCLAlign.cwl"
        },
        {
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 48000
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_sample_tag.py"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--cell-biop-summary-list"
                    },
                    "id": "#SampleTag.cwl/Cell_Biop_Summary_List"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--full-gene-list-fp"
                    },
                    "id": "#SampleTag.cwl/Full_Genes"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--feature-name-json-fp"
                    },
                    "id": "#SampleTag.cwl/Gene_Name_ID_Json"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--putative-cells-json-fp"
                    },
                    "id": "#SampleTag.cwl/Putative_Cells"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "id": "#SampleTag.cwl/Run_Metadata"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#SampleTag.cwl/num_bioproducts"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#SampleTag.cwl/num_cell_estimate"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "SampleTag/*csv"
                    },
                    "id": "#SampleTag.cwl/SampleTag_CSVs"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "SampleTag/*_Sample_Tag_Calls.csv"
                    },
                    "id": "#SampleTag.cwl/SampleTag_Calls"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "sample-tag-metrics-files.tar.gz"
                    },
                    "id": "#SampleTag.cwl/SampleTag_Metrics_Tar"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "SampleTagArchives/*zip"
                    },
                    "id": "#SampleTag.cwl/SampleTag_perTagZips"
                }
            ],
            "id": "#SampleTag.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Analyze_Reads_IG.cwl/Maximum_Threads"
                },
                {
                    "type": "long",
                    "id": "#VDJ_Analyze_Reads_IG.cwl/Num_Valid_Reads_IG"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Analyze_Reads_IG.cwl/Valid_Reads_Fastq_IG"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#VDJ_Analyze_Reads_IG.cwl/VDJ_GatherIGCalls/gatheredCalls",
                    "id": "#VDJ_Analyze_Reads_IG.cwl/gatheredCalls"
                }
            ],
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "steps": [
                {
                    "run": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/RSEC_Reads_Fastq",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/RSEC_Reads_Fastq"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Version",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/VDJ_Version"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/num_cores",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/num_cores"
                        }
                    ],
                    "out": [
                        "#VDJ_Analyze_Reads_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/igCalls"
                    ],
                    "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG"
                },
                {
                    "run": "#VDJ_GatherCalls.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Version",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_GatherIGCalls/VDJ_Version"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/igCalls",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_GatherIGCalls/theCalls"
                        }
                    ],
                    "out": [
                        "#VDJ_Analyze_Reads_IG.cwl/VDJ_GatherIGCalls/gatheredCalls"
                    ],
                    "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_GatherIGCalls"
                },
                {
                    "run": "#VDJ_Preprocess_Reads.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Analyze_Reads_IG.cwl/Maximum_Threads",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/Maximum_Threads"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_IG.cwl/Valid_Reads_Fastq_IG",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/Valid_Reads_Fastq"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_IG.cwl/Num_Valid_Reads_IG",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/num_valid_reads"
                        },
                        {
                            "valueFrom": "BCR",
                            "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/vdj_type"
                        }
                    ],
                    "out": [
                        "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/RSEC_Reads_Fastq",
                        "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/num_splits",
                        "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG/num_cores"
                    ],
                    "id": "#VDJ_Analyze_Reads_IG.cwl/VDJ_Preprocess_Reads_IG"
                }
            ],
            "id": "#VDJ_Analyze_Reads_IG.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Analyze_Reads_TCR.cwl/Maximum_Threads"
                },
                {
                    "type": "long",
                    "id": "#VDJ_Analyze_Reads_TCR.cwl/Num_Valid_Reads_TCR"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Analyze_Reads_TCR.cwl/Valid_Reads_Fastq_TCR"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_GatherTCRCalls/gatheredCalls",
                    "id": "#VDJ_Analyze_Reads_TCR.cwl/gatheredCalls"
                }
            ],
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "steps": [
                {
                    "run": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/RSEC_Reads_Fastq",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/RSEC_Reads_Fastq"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Version",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/VDJ_Version"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/num_cores",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/num_cores"
                        }
                    ],
                    "out": [
                        "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/tcrCalls"
                    ],
                    "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR"
                },
                {
                    "run": "#VDJ_GatherCalls.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Version",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_GatherTCRCalls/VDJ_Version"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/tcrCalls",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_GatherTCRCalls/theCalls"
                        }
                    ],
                    "out": [
                        "#VDJ_Analyze_Reads_TCR.cwl/VDJ_GatherTCRCalls/gatheredCalls"
                    ],
                    "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_GatherTCRCalls"
                },
                {
                    "run": "#VDJ_Preprocess_Reads.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Analyze_Reads_TCR.cwl/Maximum_Threads",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/Maximum_Threads"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_TCR.cwl/Valid_Reads_Fastq_TCR",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/Valid_Reads_Fastq"
                        },
                        {
                            "source": "#VDJ_Analyze_Reads_TCR.cwl/Num_Valid_Reads_TCR",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/num_valid_reads"
                        },
                        {
                            "valueFrom": "TCR",
                            "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/vdj_type"
                        }
                    ],
                    "out": [
                        "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/RSEC_Reads_Fastq",
                        "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/num_splits",
                        "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR/num_cores"
                    ],
                    "id": "#VDJ_Analyze_Reads_TCR.cwl/VDJ_Preprocess_Reads_TCR"
                }
            ],
            "id": "#VDJ_Analyze_Reads_TCR.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                },
                {
                    "class": "ResourceRequirement",
                    "coresMin": 1,
                    "ramMin": 1024
                }
            ],
            "baseCommand": [
                "AssembleAndAnnotate.sh"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/RSEC_Reads_Fastq"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/Read_Limit"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 3
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_pruned.csv.gz"
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/IgBlastCall"
                }
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/RSEC_Reads_Fastq"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/num_cores"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/IgBlastCall",
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/igCalls"
                }
            ],
            "steps": [
                {
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG",
                    "run": "#VDJ_Assemble_and_Annotate_Contigs.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/RSEC_Reads_Fastq",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/RSEC_Reads_Fastq"
                        },
                        {
                            "valueFrom": "75000",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/Read_Limit"
                        },
                        {
                            "source": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Version",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/IgBlastCall"
                    ],
                    "scatter": [
                        "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/RSEC_Reads_Fastq"
                    ]
                }
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/RSEC_Reads_Fastq"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/num_cores"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/IgBlastCall",
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/tcrCalls"
                }
            ],
            "steps": [
                {
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR",
                    "run": "#VDJ_Assemble_and_Annotate_Contigs.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/RSEC_Reads_Fastq",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/RSEC_Reads_Fastq"
                        },
                        {
                            "valueFrom": "75000",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/Read_Limit"
                        },
                        {
                            "source": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Version",
                            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/VDJ_Version"
                        }
                    ],
                    "out": [
                        "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/IgBlastCall"
                    ],
                    "scatter": [
                        "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/RSEC_Reads_Fastq"
                    ]
                }
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 32000
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "bdgenomics/rhapsody:4.0"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "baseCommand": [
                "mist_vdj_compile_results.py"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--seq-metrics",
                        "position": 10
                    },
                    "id": "#VDJ_Compile_Results.cwl/Seq_Metrics"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--vdj-version",
                        "position": 2
                    },
                    "id": "#VDJ_Compile_Results.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "itemSeparator": ",",
                        "prefix": "--cell-type-mapping-fp"
                    },
                    "id": "#VDJ_Compile_Results.cwl/cellTypeMapping"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--ignore",
                        "itemSeparator": ","
                    },
                    "id": "#VDJ_Compile_Results.cwl/chainsToIgnore"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 8,
                        "prefix": "--e-value-for-j"
                    },
                    "id": "#VDJ_Compile_Results.cwl/evalueJgene"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 7,
                        "prefix": "--e-value-for-v"
                    },
                    "id": "#VDJ_Compile_Results.cwl/evalueVgene"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 5
                    },
                    "id": "#VDJ_Compile_Results.cwl/igCalls"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 9,
                        "prefix": "--metadata-fp"
                    },
                    "id": "#VDJ_Compile_Results.cwl/metadata"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--putative-cells-json-fp",
                        "position": 3
                    },
                    "id": "#VDJ_Compile_Results.cwl/putativeCells"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 6
                    },
                    "id": "#VDJ_Compile_Results.cwl/tcrCalls"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "cell_order.json"
                    },
                    "id": "#VDJ_Compile_Results.cwl/Cell_Order"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "scirpy compatible h5ad output that contains all the contig information present in the AIRR file",
                    "outputBinding": {
                        "glob": "*_VDJ.h5ad"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjAnnData"
                },
                {
                    "doc": "VDJ data per cell, with distribution based error correction",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_perCell.csv"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjCellsDatatable"
                },
                {
                    "doc": "VDJ data per cell, including non-putative cells, no error correction applied",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_perCell_uncorrected.csv.gz"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjCellsDatatableUncorrected"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_DBEC_images.tar.gz"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjDbecFilterImages"
                },
                {
                    "doc": "AIRR compatible output that only reports the Dominant contigs, counts are DBEC corrected",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_Dominant_Contigs_AIRR.tsv"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjDominantContigsAIRR"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_VDJ_metrics.csv"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjMetricsCsv"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "vdj-metrics-files.tar.gz"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjMetricsTar"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "AIRR compatible output that reports all the contigs, counts are not DBEC corrected",
                    "outputBinding": {
                        "glob": "*_VDJ_Unfiltered_Contigs_AIRR.tsv.gz"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjUnfilteredContigsAIRR"
                }
            ],
            "id": "#VDJ_Compile_Results.cwl"
        },
        {
            "doc": "VDJ_GatherCalls collect the outputs from the multi-processed VDJ step into one file.\n",
            "class": "Workflow",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_GatherCalls.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#VDJ_GatherCalls.cwl/theCalls"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gatheredCalls",
                    "id": "#VDJ_GatherCalls.cwl/gatheredCalls"
                }
            ],
            "steps": [
                {
                    "in": [
                        {
                            "source": "#VDJ_GatherCalls.cwl/theCalls",
                            "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/theCalls"
                        }
                    ],
                    "out": [
                        "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gatheredCalls"
                    ],
                    "run": {
                        "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/run/gather_IgBlast",
                        "cwlVersion": "v1.2",
                        "class": "CommandLineTool",
                        "hints": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:4.0",
                                "class": "DockerRequirement"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            },
                            {
                                "class": "ShellCommandRequirement"
                            }
                        ],
                        "inputs": [
                            {
                                "type": [
                                    "null",
                                    {
                                        "type": "array",
                                        "items": "File"
                                    }
                                ],
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/run/gather_IgBlast/theCalls"
                            }
                        ],
                        "arguments": [
                            {
                                "shellQuote": false,
                                "valueFrom": "${\n  if (!inputs.theCalls[0] ) {\n    return (\"echo \\\"No outputs from IgBlast detected in VDJ_GatherCalls\\\"\")\n  }\n  var inputFiles = \"\"\n  if (!inputs.theCalls[0].path.split(\"_PruneIgBlast\")[1]){\n    inputFiles = \"zcat\"\n    for (var i = 0; i < inputs.theCalls.length; i++) {\n      inputFiles += \" \" + inputs.theCalls[i].path\n    }\n    inputFiles += \" | \"\n  } else {\n    inputFiles = \"zcat \" + inputs.theCalls[0].path.split(\"VDJ\")[0] + \"*\" + inputs.theCalls[0].path.split(\"_PruneIgBlast\")[1].split(\"_Number_\")[0] + \"_Number_*.csv.gz | \"\n  }\n  var outputFileName = \"\\\"gzip > \" + inputs.theCalls[0].nameroot.split(\"_Number_\")[0] + \"_constant_region_called_pruned.csv.gz\" + \"\\\"\"\n  var awkCommand =  \"awk \\'NR==1{F=$1;print | \" + outputFileName + \" } $1!=F { print | \" + outputFileName + \" }\\' \"\n  var outputCommand = inputFiles + awkCommand\n  return (outputCommand)\n}"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "outputBinding": {
                                    "glob": "*_constant_region_called_pruned.csv.gz",
                                    "outputEval": "${\n  if (self.size == 0) {\n    throw(\"No outputs from IgBlast detected in VDJ_GatherCalls!\");\n  } else {\n    return(self);\n  }\n}"
                                },
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/run/gather_IgBlast/gatheredCalls"
                            }
                        ]
                    },
                    "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls"
                }
            ],
            "id": "#VDJ_GatherCalls.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/Maximum_Threads"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/Valid_Reads_Fastq"
                },
                {
                    "type": "long",
                    "id": "#VDJ_Preprocess_Reads.cwl/num_valid_reads"
                },
                {
                    "type": "string",
                    "id": "#VDJ_Preprocess_Reads.cwl/vdj_type"
                }
            ],
            "outputs": [
                {
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/RSEC_Reads_Fastq",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#VDJ_Preprocess_Reads.cwl/RSEC_Reads_Fastq"
                },
                {
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_cores",
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/num_cores"
                },
                {
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_splits",
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/num_splits"
                }
            ],
            "requirements": [
                {
                    "envDef": [
                        {
                            "envValue": "8",
                            "envName": "CORES_ALLOCATED_PER_CWL_PROCESS"
                        }
                    ],
                    "class": "EnvVarRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "steps": [
                {
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads",
                    "hints": [
                        {
                            "class": "ResourceRequirement",
                            "coresMin": 8
                        }
                    ],
                    "run": "#VDJ_RSEC_Reads.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Valid_Reads",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/Valid_Reads"
                        },
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_splits",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/num_splits"
                        },
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/num_valid_reads",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/num_valid_reads"
                        }
                    ],
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/RSEC_Reads_Fastq"
                    ]
                },
                {
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads",
                    "hints": [
                        {
                            "class": "ResourceRequirement",
                            "coresMin": 8
                        }
                    ],
                    "run": "#VDJ_Trim_Reads.cwl",
                    "in": [
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/Valid_Reads_Fastq",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Valid_Reads_Fastq"
                        }
                    ],
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Valid_Reads",
                        "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Trim_Report"
                    ]
                },
                {
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits",
                    "in": [
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/Maximum_Threads",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/Maximum_Threads"
                        },
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/num_valid_reads",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_valid_reads"
                        },
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/vdj_type",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/vdj_type"
                        }
                    ],
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_splits",
                        "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_cores"
                    ],
                    "run": {
                        "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits",
                        "cwlVersion": "v1.2",
                        "class": "ExpressionTool",
                        "inputs": [
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/Maximum_Threads"
                            },
                            {
                                "type": [
                                    "null",
                                    "long"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/num_valid_reads"
                            },
                            {
                                "type": "string",
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/vdj_type"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/num_cores"
                            },
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/run/determine_num_splits/num_splits"
                            }
                        ],
                        "expression": "${\n  var num_splits = 64;\n  var max_threads = parseInt(inputs.Maximum_Threads);\n  if (!isNaN(max_threads)) {\n    num_splits = parseInt(Math.max(max_threads, 8) * 0.7);\n  }\n  return ({\"num_splits\": num_splits, \"num_cores\": num_splits});\n}"
                    }
                }
            ],
            "id": "#VDJ_Preprocess_Reads.cwl"
        },
        {
            "hints": [
                {
                    "dockerPull": "bdgenomics/rhapsody:4.0",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": "mist_vdj_rsec_reads.py",
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_RSEC_Reads.cwl/VDJ_Version"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "prefix": "--vdj-valid-reads",
                        "itemSeparator": ","
                    },
                    "id": "#VDJ_RSEC_Reads.cwl/Valid_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--num-splits"
                    },
                    "id": "#VDJ_RSEC_Reads.cwl/num_splits"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#VDJ_RSEC_Reads.cwl/num_valid_reads"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*RSEC_Reads_Fastq_*.tar.gz"
                    },
                    "id": "#VDJ_RSEC_Reads.cwl/RSEC_Reads_Fastq"
                }
            ],
            "id": "#VDJ_RSEC_Reads.cwl"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#VDJ_Settings.cwl/_VDJ_Version"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_VGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_Version"
                }
            ],
            "expression": "${\n  var vdjVersion = null;\n  if (!inputs._VDJ_Version) {\n    vdjVersion = null;}\n  else {\n    var _VDJ_Version = inputs._VDJ_Version.toLowerCase();\n    if (_VDJ_Version === \"human\" || _VDJ_Version === \"hs\" || _VDJ_Version === \"human vdj - bcr and tcr\") {\n      vdjVersion = \"human\";\n    } else if (_VDJ_Version === \"humanbcr\" || _VDJ_Version === \"human vdj - bcr only\") {\n      vdjVersion = \"humanBCR\";\n    } else if (_VDJ_Version === \"humantcr\" || _VDJ_Version === \"human vdj - tcr only\") {\n      vdjVersion = \"humanTCR\";\n    } else if (_VDJ_Version === \"mouse\" || _VDJ_Version === \"mm\" || _VDJ_Version === \"mouse vdj - bcr and tcr\") {\n      vdjVersion = \"mouse\";\n    } else if (_VDJ_Version === \"mousebcr\" || _VDJ_Version === \"mouse vdj - bcr only\") {\n      vdjVersion = \"mouseBCR\";\n    } else if (_VDJ_Version === \"mousetcr\" || _VDJ_Version === \"mouse vdj - tcr only\") {\n      vdjVersion = \"mouseTCR\";\n    } else {\n      vdjVersion = inputs._VDJ_Version;\n    }\n  }\n\n  return ({\n  VDJ_Version: vdjVersion,\n  })\n}",
            "id": "#VDJ_Settings.cwl"
        },
        {
            "hints": [
                {
                    "dockerPull": "bdgenomics/rhapsody:4.0",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": "VDJ_Trim_Reads.sh",
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Trim_Reads.cwl/VDJ_Version"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#VDJ_Trim_Reads.cwl/Valid_Reads_Fastq"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "cutadapt.log"
                    },
                    "id": "#VDJ_Trim_Reads.cwl/Trim_Report"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "*vdjtxt.gz"
                    },
                    "id": "#VDJ_Trim_Reads.cwl/Valid_Reads"
                }
            ],
            "id": "#VDJ_Trim_Reads.cwl"
        },
        {
            "hints": [
                {
                    "dockerPull": "bdgenomics/rhapsody:4.0",
                    "class": "DockerRequirement"
                }
            ],
            "class": "CommandLineTool",
            "baseCommand": [
                "mist_print_version.py"
            ],
            "stdout": "output.txt",
            "inputs": [],
            "outputs": [
                {
                    "type": "string",
                    "outputBinding": {
                        "glob": "output.txt",
                        "loadContents": true,
                        "outputEval": "$(self[0].contents)"
                    },
                    "id": "#Version.cwl/version"
                }
            ],
            "id": "#Version.cwl"
        }
    ],
    "cwlVersion": "v1.2"
}