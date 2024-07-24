#!/usr/bin/env cwl-runner
{
    "cwlVersion": "v1.0",
    "$graph": [
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--annot-r1",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AddtoBam.cwl/Annotation_R1"
                },
                {
                    "inputBinding": {
                        "prefix": "--cell-order"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/Cell_Order"
                },
                {
                    "inputBinding": {
                        "prefix": "--annot-mol-file"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/Molecular_Annotation"
                },
                {
                    "inputBinding": {
                        "prefix": "--r2-bam"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/R2_Bam"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--tag-calls"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AddtoBam.cwl/Tag_Calls"
                },
                {
                    "inputBinding": {
                        "prefix": "--target-gene-mapping"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AddtoBam.cwl/Target_Gene_Mapping"
                }
            ],
            "requirements": [
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "Annotated_mapping_R2.BAM"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/Annotated_Bam"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_add_to_bam.py"
            ],
            "class": "CommandLineTool",
            "id": "#AddtoBam.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignR2.cwl/Extra_Seqs"
                },
                {
                    "inputBinding": {
                        "prefix": "--index"
                    },
                    "type": "File",
                    "id": "#AlignR2.cwl/Index"
                },
                {
                    "inputBinding": {
                        "prefix": "--r2-fastqs",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AlignR2.cwl/R2"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#AlignR2.cwl/Run_Metadata"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "envDef": [
                        {
                            "envName": "CORES_ALLOCATED_PER_CWL_PROCESS",
                            "envValue": "$(String(runtime.cores))"
                        }
                    ],
                    "class": "EnvVarRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*zip"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AlignR2.cwl/Alignments"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AlignR2.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_align_R2.py"
            ],
            "class": "CommandLineTool",
            "id": "#AlignR2.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--umi-option"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateMolecules.cwl/AbSeq_UMI"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--use-dbec"
                    },
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#AnnotateMolecules.cwl/Use_DBEC"
                },
                {
                    "inputBinding": {
                        "prefix": "--valid-annot"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/Valids"
                }
            ],
            "requirements": [
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_GeneStatus.csv.*"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/Gene_Status_List"
                },
                {
                    "outputBinding": {
                        "glob": "stats.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).max_count)\n"
                    },
                    "type": "int",
                    "id": "#AnnotateMolecules.cwl/Max_Count"
                },
                {
                    "outputBinding": {
                        "glob": "*_Annotation_Molecule.csv.*"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/Mol_Annot_List"
                },
                {
                    "outputBinding": {
                        "glob": "stats.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).total_molecules)\n"
                    },
                    "type": "int",
                    "id": "#AnnotateMolecules.cwl/Total_Molecules"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_annotate_molecules.py"
            ],
            "class": "CommandLineTool",
            "id": "#AnnotateMolecules.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--filter-metrics",
                        "itemSeparator": ","
                    },
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#AnnotateR1.cwl/Filter_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--R1"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/R1"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/Run_Metadata"
                }
            ],
            "requirements": [
                {
                    "ramMin": 2000,
                    "class": "ResourceRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_Annotation_R1.csv.gz"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/Annotation_R1"
                },
                {
                    "outputBinding": {
                        "glob": "*_R1_error_count_table.npy"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/R1_error_count_table"
                },
                {
                    "outputBinding": {
                        "glob": "*_R1_read_count_breakdown.json"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/R1_read_count_breakdown"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_annotate_R1.py"
            ],
            "class": "CommandLineTool",
            "id": "#AnnotateR1.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateR2.cwl/Extra_Seqs"
                },
                {
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateR2.cwl/GTF_Annotation"
                },
                {
                    "inputBinding": {
                        "prefix": "--R2-zip"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/R2_zip"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--transcript-length"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateR2.cwl/Transcript_Length"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*Annotation_R2.csv.gz"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/Annot_R2"
                },
                {
                    "outputBinding": {
                        "glob": "*-annot.gtf"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateR2.cwl/GTF"
                },
                {
                    "outputBinding": {
                        "glob": "*mapping_R2.BAM"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/R2_Bam"
                },
                {
                    "outputBinding": {
                        "glob": "*_picard_quality_metrics.csv.gz"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/R2_Quality_Metrics"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_annotate_R2.py"
            ],
            "class": "CommandLineTool",
            "id": "#AnnotateR2.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--umi-option"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateReads.cwl/AbSeq_UMI"
                },
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/Extra_Seqs"
                },
                {
                    "type": {
                        "items": [
                            "null",
                            "File"
                        ],
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/Filter_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateReads.cwl/Putative_Cell_Call"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/R1_Annotation"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/R1_error_count_table"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/R1_read_count_breakdown"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/R2_Annotation"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/R2_Quality_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#AnnotateReads.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--target-gene-mapping"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/Target_Gene_Mapping"
                }
            ],
            "requirements": [
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "writable": false,
                            "entry": "${\n    function getPaths(inputs, attribute) {\n      var fp_arr = []\n      for (var i = 0; i < inputs[attribute].length; i++)\n      {\n          fp_arr.push(inputs[attribute][i].path);\n      }\n      return fp_arr;\n    }\n    var paths = {}\n    paths['annotR1'] = getPaths(inputs, 'R1_Annotation')\n    paths['R1_error_count_table'] = getPaths(inputs, 'R1_error_count_table')\n    paths['R1_read_count_breakdown'] = getPaths(inputs, 'R1_read_count_breakdown')\n    paths['annotR2'] = getPaths(inputs, 'R2_Annotation')\n    paths['r2_quality_metrics_fps'] = getPaths(inputs, 'R2_Quality_Metrics')\n    if(inputs.Filter_Metrics[0] != null){\n        paths['filtering_stat_files'] = getPaths(inputs, 'Filter_Metrics')\n    }\n    var paths_json = JSON.stringify(paths);\n    return paths_json;\n}",
                            "entryname": "manifest.json"
                        }
                    ]
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "envDef": [
                        {
                            "envName": "CORES_ALLOCATED_PER_CWL_PROCESS",
                            "envValue": "4"
                        }
                    ],
                    "class": "EnvVarRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_Annotation_Read.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/Annotation_Read"
                },
                {
                    "outputBinding": {
                        "glob": "*read1_error_rate_archive*"
                    },
                    "type": "File",
                    "id": "#AnnotateReads.cwl/Read1_error_rate"
                },
                {
                    "outputBinding": {
                        "glob": "*_SeqMetrics.csv.gz"
                    },
                    "type": "File",
                    "id": "#AnnotateReads.cwl/Seq_Metrics"
                },
                {
                    "outputBinding": {
                        "glob": "*Sorted_Valid_Reads.csv.*"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/Valid_Reads"
                },
                {
                    "outputBinding": {
                        "glob": "num_vdj_reads.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 0; } return parseInt(JSON.parse(self[0].contents).BCR); }"
                    },
                    "type": "int",
                    "id": "#AnnotateReads.cwl/num_valid_ig_reads"
                },
                {
                    "outputBinding": {
                        "glob": "num_vdj_reads.json",
                        "loadContents": true,
                        "outputEval": "${ if (!self[0]) { return 0; } return parseInt(JSON.parse(self[0].contents).TCR); }"
                    },
                    "type": "int",
                    "id": "#AnnotateReads.cwl/num_valid_tcr_reads"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AnnotateReads.cwl/output"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_IG_Valid_Reads.fastq.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/validIgReads"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_TCR_Valid_Reads.fastq.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/validTcrReads"
                }
            ],
            "baseCommand": [
                "mist_annotate_reads.py"
            ],
            "class": "CommandLineTool",
            "id": "#AnnotateReads.cwl"
        },
        {
            "inputs": [
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#BundleLogs.cwl/log_files"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "id": "#BundleLogs.cwl/logs_dir"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */\n  function uuid() {\n    var uuid = \"\", i, random;\n    for (i = 0; i < 32; i++) {\n      random = Math.random() * 16 | 0;\n      if (i == 8 || i == 12 || i == 16 || i == 20) {\n        uuid += \"-\";\n      }\n      uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);\n    }\n    return uuid;\n  }\n  var listing = [];\n  for (var i = 0; i < inputs.log_files.length; i++) {\n    var log_file = inputs.log_files[i];\n    log_file.basename = uuid() + \"-\" + log_file.basename;\n    listing.push(log_file);\n  }\n  return ({\n    logs_dir: {\n      class: \"Directory\",\n      basename: \"Logs\",\n      listing: listing\n    }\n  });\n}",
            "id": "#BundleLogs.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 0
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Cell_Classifier.cwl/molsPerCellMatrix"
                }
            ],
            "requirements": [
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*cell_type_experimental.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Cell_Classifier.cwl/cellTypePredictions"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#Cell_Classifier.cwl/log"
                }
            ],
            "baseCommand": [
                "mist_cell_classifier.py"
            ],
            "class": "CommandLineTool",
            "id": "#Cell_Classifier.cwl"
        },
        {
            "inputs": [
                {
                    "doc": "The minimum size (megabytes) of a file that should get split into chunks of a size designated in NumRecordsPerSplit\n",
                    "inputBinding": {
                        "prefix": "--min-split-size"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CheckFastqs.cwl/MinChunkSize"
                },
                {
                    "inputBinding": {
                        "prefix": "--reads",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#CheckFastqs.cwl/Reads"
                },
                {
                    "inputBinding": {
                        "prefix": "--subsample"
                    },
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#CheckFastqs.cwl/Subsample"
                },
                {
                    "inputBinding": {
                        "prefix": "--subsample-seed"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CheckFastqs.cwl/Subsample_Seed"
                },
                {
                    "inputBinding": {
                        "prefix": "--subsample-seed"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CheckFastqs.cwl/UserInputSubsampleSeed"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "CheckFastqs does several quality control routines including: (1) ensuring that read pair file names are formatted correctly and contain a read pair mate; (2) disambiguating the \"Subsample Reads\" input and; (3) if not provided, generating a subsampling seed that the downstream instances can use.\n",
            "baseCommand": [
                "mist_check_fastqs.py"
            ],
            "id": "#CheckFastqs.cwl",
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "bead_version.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).BeadVersion)\n"
                    },
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/Bead_Version/Library"
                                },
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/Bead_Version/bead_version"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#CheckFastqs.cwl/Bead_Version"
                },
                {
                    "outputBinding": {
                        "glob": "fastq_read_pairs.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).fastq_read_pairs)\n"
                    },
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/filename"
                                },
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/readFlag"
                                },
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/readPairId"
                                },
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/library"
                                },
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/beadVersion"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#CheckFastqs.cwl/FastqReadPairs"
                },
                {
                    "outputBinding": {
                        "glob": "files_to_skip_split_and_subsample.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).files_to_skip_split_and_subsample)\n"
                    },
                    "type": {
                        "items": "string",
                        "type": "array"
                    },
                    "id": "#CheckFastqs.cwl/FilesToSkipSplitAndSubsample"
                },
                {
                    "outputBinding": {
                        "glob": "fastq_read_pairs.json",
                        "loadContents": true,
                        "outputEval": "${\n  var obj = JSON.parse(self[0].contents);\n  var libraries = [];\n  var pairs = obj.fastq_read_pairs\n  for (var i in pairs){\n    if (pairs[i][\"readFlag\"] == \"R1\"){\n      if (libraries.indexOf(pairs[i][\"library\"]) == -1){ \n        libraries.push(pairs[i][\"library\"]);\n      }\n    }\n  }\n  libraries.sort();\n  return(libraries.toString())\n}\n"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#CheckFastqs.cwl/Libraries"
                },
                {
                    "outputBinding": {
                        "outputEval": "${  \n  var reads = []; \n  var files = inputs.Reads\n  for (var i in files){\n      reads.push(files[i][\"basename\"]);\n  }\n  reads.sort();\n  return(reads)\n}\n"
                    },
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#CheckFastqs.cwl/ReadsList"
                },
                {
                    "outputBinding": {
                        "glob": "subsampling_info.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).subsampling_seed)\n"
                    },
                    "type": "int",
                    "id": "#CheckFastqs.cwl/SubsampleSeed"
                },
                {
                    "outputBinding": {
                        "glob": "subsampling_info.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).subsampling_ratio)\n"
                    },
                    "type": "float",
                    "id": "#CheckFastqs.cwl/SubsamplingRatio"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#CheckFastqs.cwl/log"
                }
            ],
            "class": "CommandLineTool"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--abseq-reference",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#CheckReference.cwl/AbSeq_Reference"
                },
                {
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CheckReference.cwl/Putative_Cell_Call"
                },
                {
                    "inputBinding": {
                        "prefix": "--reference",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#CheckReference.cwl/Reference"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#CheckReference.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--supplemental-reference",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#CheckReference.cwl/Supplemental_Reference"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "combined_extra_seq.fasta"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Extra_Seqs"
                },
                {
                    "outputBinding": {
                        "glob": "full-gene-list.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Full_Genes"
                },
                {
                    "outputBinding": {
                        "glob": "*gtf",
                        "outputEval": "${\n    // get the WTA modified GTF with extra seqs\n    if (self.length == 1) {\n        return self;\n    // there is no modified GTF\n    } else if (self.length == 0) {\n        // if Reference is null (i.e. AbSeq_Reference only), return no GTF\n        if (inputs.Reference === null) {\n            return null;\n        } else {\n            // get the original WTA GTF without extra seqs\n            for (var i = 0; i < inputs.Reference.length; i++) {\n                if (inputs.Reference[i].basename.toLowerCase().indexOf('gtf') !== -1) {\n                    return inputs.Reference[i];\n                }\n            }\n            // return no GTF for Targeted\n            return null\n        }\n    }\n}\n"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/GTF"
                },
                {
                    "outputBinding": {
                        "glob": "*-annot.*",
                        "outputEval": "${\n    if (self.length == 1) { // Targeted\n        return self;\n    } else if (self.length == 0){ // WTA without extra seqs or targets\n        for (var i = 0; i < inputs.Reference.length; i++) {\n            if (inputs.Reference[i].basename.toLowerCase().indexOf('tar.gz') !== -1) {\n                return inputs.Reference[i];\n            }\n        }\n        return null\n    }\n}\n"
                    },
                    "type": "File",
                    "id": "#CheckReference.cwl/Index"
                },
                {
                    "outputBinding": {
                        "glob": "target-gene.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Target_Gene_Mapping"
                },
                {
                    "outputBinding": {
                        "glob": "transcript_length.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Transcript_Length"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#CheckReference.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_check_references.py"
            ],
            "class": "CommandLineTool",
            "id": "#CheckReference.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--cell-order"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/Cell_Order"
                },
                {
                    "inputBinding": {
                        "prefix": "--dense-data-table"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#DensetoSparse.cwl/Dense_Data_Table"
                },
                {
                    "inputBinding": {
                        "prefix": "--gene-list"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/Gene_List"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/Run_Metadata"
                }
            ],
            "requirements": [
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*.csv.gz"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/Data_Tables"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_dense_to_sparse.py"
            ],
            "class": "CommandLineTool",
            "id": "#DensetoSparse.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#DensetoSparseFile.cwl/GDT_cell_order"
                }
            ],
            "requirements": [
            ],
            "stdout": "cell_order.json",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#DensetoSparseFile.cwl/Cell_Order"
                }
            ],
            "baseCommand": "cat",
            "id": "#DensetoSparseFile.cwl",
            "class": "CommandLineTool"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--full-gene-list"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Full_Genes"
                },
                {
                    "inputBinding": {
                        "prefix": "--gene-status",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Gene_Status_List"
                },
                {
                    "inputBinding": {
                        "prefix": "--max-count",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "int",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Max_Count"
                },
                {
                    "inputBinding": {
                        "prefix": "--mol-annot",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Molecule_Annotation_List"
                },
                {
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#GetDataTable.cwl/Putative_Cell_Call"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--seq-metrics"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Seq_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--tag-names",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#GetDataTable.cwl/Tag_Names"
                },
                {
                    "type": {
                        "items": "int",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Total_Molecules"
                }
            ],
            "requirements": [
                {
                    "ramMin": "${return Math.min(Math.max(parseInt(inputs.Total_Molecules.reduce(function(a, b) { return a + b; }, 0) / 4000), 32000), 768000);}",
                    "class": "ResourceRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "metrics-files.tar.gz"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Annot_Files"
                },
                {
                    "outputBinding": {
                        "glob": "Annotations/*_Bioproduct_Stats.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Bioproduct_Stats"
                },
                {
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*.png"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#GetDataTable.cwl/Cell_Label_Filter"
                },
                {
                    "outputBinding": {
                        "glob": "cell_order.json"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Cell_Order"
                },
                {
                    "outputBinding": {
                        "glob": "*_Annotation_Molecule_corrected.csv.gz"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Corrected_Molecular_Annotation"
                },
                {
                    "outputBinding": {
                        "glob": "*PerCell_Dense.csv.gz"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Dense_Data_Tables"
                },
                {
                    "outputBinding": {
                        "glob": "*PerCell_Unfiltered_Dense.csv.gz"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Dense_Data_Tables_Unfiltered"
                },
                {
                    "outputBinding": {
                        "glob": "*_Expression_Data.st.gz"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Expression_Data"
                },
                {
                    "outputBinding": {
                        "glob": "*_Expression_Data_Unfiltered.st.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Expression_Data_Unfiltered"
                },
                {
                    "outputBinding": {
                        "glob": "gene_list.json"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Gene_List"
                },
                {
                    "outputBinding": {
                        "glob": "Annotations/*_Annotation_Molecule.csv.gz"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Molecular_Annotation"
                },
                {
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*_Protein_Aggregates_Experimental.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Protein_Aggregates_Experimental"
                },
                {
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*_Putative_Cells_Origin.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Putative_Cells_Origin"
                },
                {
                    "outputBinding": {
                        "glob": "Annotations/*_Annotation_Molecule_Trueno.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Tag_Annotation"
                },
                {
                    "outputBinding": {
                        "glob": "Trueno/*_Calls.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Tag_Calls"
                },
                {
                    "outputBinding": {
                        "glob": "Trueno/*csv"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#GetDataTable.cwl/Trueno_out"
                },
                {
                    "outputBinding": {
                        "glob": "Trueno/*zip"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#GetDataTable.cwl/Trueno_zip"
                },
                {
                    "outputBinding": {
                        "glob": "Annotations/*_UMI_Adjusted_CellLabel_Stats.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/UMI_Adjusted_CellLabel_Stats"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_get_datatables.py"
            ],
            "class": "CommandLineTool",
            "id": "#GetDataTable.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    },
                    "type": "File",
                    "id": "#IndexBAM.cwl/BamFile"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "samtools_index.log",
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*.bai"
                    },
                    "type": "File",
                    "id": "#IndexBAM.cwl/Index"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#IndexBAM.cwl/log"
                }
            ],
            "baseCommand": [
                "samtools",
                "index"
            ],
            "id": "#IndexBAM.cwl",
            "arguments": [
                {
                    "position": 2,
                    "valueFrom": "${\n    return inputs.BamFile.basename + \".bai\"\n}"
                }
            ],
            "class": "CommandLineTool"
        },
        {
            "inputs": [],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
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
                    "id": "#InternalSettings.cwl/Barcode_Num"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#InternalSettings.cwl/Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/Label_Version"
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
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Read_Filter_Off"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#InternalSettings.cwl/Seq_Run"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/Subsample_Tags"
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
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Use_DBEC"
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
            "class": "ExpressionTool",
            "expression": "${\n  var internalInputs = [\n    '_Label_Version',\n    '_Read_Filter_Off',\n    '_Barcode_Num',\n    '_Seq_Run',\n    '_AbSeq_UMI',\n    '_Use_DBEC',\n    '_Extra_Seqs',\n    '_MinChunkSize',\n    '_NumRecordsPerSplit',\n    '_Target_analysis',\n    '_Subsample_Tags',\n    '_VDJ_VGene_Evalue',\n    '_VDJ_JGene_Evalue',\n  ];\n  var internalOutputs = {}\n  for (var i = 0; i < internalInputs.length; i++) {\n    var internalInput = internalInputs[i];\n    var internalOutput = internalInput.slice(1); // remove leading underscore\n    if (inputs.hasOwnProperty(internalInput)) {\n      internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output\n    } else {\n      internalOutputs[internalOutput] = null; // if input not specified, provide a null\n    }\n  }\n  return internalOutputs;\n}",
            "id": "#InternalSettings.cwl"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/AbSeq_Reference",
                    "label": "AbSeq Reference"
                },
                {
                    "doc": "Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set.",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#main/Basic_Algo_Only",
                    "label": "Disable Refined Putative Cell Calling"
                },
                {
                    "doc": "Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count",
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/Exact_Cell_Count",
                    "label": "Exact Cell Count"
                },
                {
                    "doc": "Specify the data to be used for putative cell calling. mRNA is the default selected option. AbSeq (Experimental) is for troubleshooting only.",
                    "type": [
                        "null",
                        {
                            "symbols": [
                                "#main/Putative_Cell_Call/Putative_Cell_Call/mRNA",
                                "#main/Putative_Cell_Call/Putative_Cell_Call/AbSeq_Experimental"
                            ],
                            "type": "enum",
                            "name": "#main/Putative_Cell_Call/Putative_Cell_Call"
                        }
                    ],
                    "id": "#main/Putative_Cell_Call",
                    "label": "Putative Cell Calling"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#main/Reads",
                    "label": "Reads"
                },
                {
                    "doc": "A fasta file containing the mRNA panel amplicon targets used in the experiment",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Reference",
                    "label": "Reference"
                },
                {
                    "doc": "This is a name for output files, for example Experiment1_Metrics_Summary.csv. Default if left empty is to name run based on a library. Any non-alpha numeric characters will be changed to a hyphen.",
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#main/Run_Name",
                    "label": "Run Name"
                },
                {
                    "doc": "The sample multiplexing kit version.  This option should only be set for a multiplexed experiment.",
                    "type": [
                        "null",
                        {
                            "symbols": [
                                "#main/Sample_Tags_Version/Sample_Tags_Version/human",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/hs",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/mouse",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/mm",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/custom"
                            ],
                            "type": "enum",
                            "name": "#main/Sample_Tags_Version/Sample_Tags_Version"
                        }
                    ],
                    "id": "#main/Sample_Tags_Version",
                    "label": "Sample Tags Version"
                },
                {
                    "doc": "Any number of reads >1 or a fraction between 0 < n < 1 to indicate the percentage of reads to subsample.\n",
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#main/Subsample",
                    "label": "Subsample Reads"
                },
                {
                    "doc": "For use when replicating a previous subsampling run only. Obtain the seed generated from the log file for the SplitFastQ node.\n",
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/Subsample_seed",
                    "label": "Subsample Seed"
                },
                {
                    "doc": "Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Should be alpha numeric, with + - and _ allowed. Any special characters: &, (), [], {}, <>, ?, | will be corrected to underscores. \n",
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Tag_Names",
                    "label": "Tag Names"
                },
                {
                    "doc": "The VDJ species and chain types.  This option should only be set for VDJ experiment.",
                    "type": [
                        "null",
                        {
                            "symbols": [
                                "#main/VDJ_Version/VDJ_Version/human",
                                "#main/VDJ_Version/VDJ_Version/hs",
                                "#main/VDJ_Version/VDJ_Version/mouse",
                                "#main/VDJ_Version/VDJ_Version/mm",
                                "#main/VDJ_Version/VDJ_Version/humanBCR",
                                "#main/VDJ_Version/VDJ_Version/humanTCR",
                                "#main/VDJ_Version/VDJ_Version/mouseBCR",
                                "#main/VDJ_Version/VDJ_Version/mouseTCR"
                            ],
                            "type": "enum",
                            "name": "#main/VDJ_Version/VDJ_Version"
                        }
                    ],
                    "id": "#main/VDJ_Version",
                    "label": "VDJ Species Version"
                }
            ],
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
            "doc": "The BD Rhapsody\u2122 assays are used to create sequencing libraries from single cell transcriptomes.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files and a reference file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
            "label": "BD Rhapsody\u2122 Targeted Analysis Pipeline",
            "steps": [
                {
                    "run": "#AddtoBam.cwl",
                    "scatter": [
                        "#main/AddtoBam/R2_Bam"
                    ],
                    "in": [
                        {
                            "source": "#main/AnnotateR1/Annotation_R1",
                            "id": "#main/AddtoBam/Annotation_R1"
                        },
                        {
                            "source": "#main/Dense_to_Sparse_File/Cell_Order",
                            "id": "#main/AddtoBam/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Corrected_Molecular_Annotation",
                            "id": "#main/AddtoBam/Molecular_Annotation"
                        },
                        {
                            "source": "#main/AnnotateR2/R2_Bam",
                            "id": "#main/AddtoBam/R2_Bam"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AddtoBam/Run_Metadata"
                        },
                        {
                            "source": "#main/GetDataTable/Tag_Calls",
                            "id": "#main/AddtoBam/Tag_Calls"
                        },
                        {
                            "source": "#main/CheckReference/Target_Gene_Mapping",
                            "id": "#main/AddtoBam/Target_Gene_Mapping"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 16000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AddtoBam",
                    "out": [
                        "#main/AddtoBam/Annotated_Bam",
                        "#main/AddtoBam/output"
                    ]
                },
                {
                    "run": "#AlignR2.cwl",
                    "out": [
                        "#main/AlignR2/Alignments",
                        "#main/AlignR2/output"
                    ],
                    "requirements": [
                        {
                            "coresMin": 8,
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AlignR2",
                    "in": [
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AlignR2/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/Index",
                            "id": "#main/AlignR2/Index"
                        },
                        {
                            "source": "#main/QualityFilterOuter/R2",
                            "id": "#main/AlignR2/R2"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AlignR2/Run_Metadata"
                        }
                    ]
                },
                {
                    "run": "#AnnotateMolecules.cwl",
                    "scatter": [
                        "#main/AnnotateMolecules/Valids"
                    ],
                    "in": [
                        {
                            "source": "#main/Internal_Settings/AbSeq_UMI",
                            "id": "#main/AnnotateMolecules/AbSeq_UMI"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateMolecules/Run_Metadata"
                        },
                        {
                            "source": "#main/Internal_Settings/Use_DBEC",
                            "id": "#main/AnnotateMolecules/Use_DBEC"
                        },
                        {
                            "source": "#main/AnnotateReads/Valid_Reads",
                            "id": "#main/AnnotateMolecules/Valids"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 32000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AnnotateMolecules",
                    "out": [
                        "#main/AnnotateMolecules/Mol_Annot_List",
                        "#main/AnnotateMolecules/Gene_Status_List",
                        "#main/AnnotateMolecules/Max_Count",
                        "#main/AnnotateMolecules/Total_Molecules",
                        "#main/AnnotateMolecules/output"
                    ]
                },
                {
                    "id": "#main/AnnotateR1",
                    "out": [
                        "#main/AnnotateR1/Annotation_R1",
                        "#main/AnnotateR1/R1_error_count_table",
                        "#main/AnnotateR1/R1_read_count_breakdown",
                        "#main/AnnotateR1/output"
                    ],
                    "run": "#AnnotateR1.cwl",
                    "scatter": [
                        "#main/AnnotateR1/R1"
                    ],
                    "in": [
                        {
                            "source": "#main/QualityFilterOuter/Filter_Metrics",
                            "id": "#main/AnnotateR1/Filter_Metrics"
                        },
                        {
                            "source": "#main/QualityFilterOuter/R1",
                            "id": "#main/AnnotateR1/R1"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateR1/Run_Metadata"
                        }
                    ]
                },
                {
                    "run": "#AnnotateR2.cwl",
                    "scatter": [
                        "#main/AnnotateR2/R2_zip"
                    ],
                    "in": [
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AnnotateR2/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/AnnotateR2/GTF_Annotation"
                        },
                        {
                            "source": "#main/AlignR2/Alignments",
                            "id": "#main/AnnotateR2/R2_zip"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateR2/Run_Metadata"
                        },
                        {
                            "source": "#main/CheckReference/Transcript_Length",
                            "id": "#main/AnnotateR2/Transcript_Length"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AnnotateR2",
                    "out": [
                        "#main/AnnotateR2/Annot_R2",
                        "#main/AnnotateR2/R2_Bam",
                        "#main/AnnotateR2/GTF",
                        "#main/AnnotateR2/output",
                        "#main/AnnotateR2/R2_Quality_Metrics"
                    ]
                },
                {
                    "run": "#AnnotateReads.cwl",
                    "out": [
                        "#main/AnnotateReads/Seq_Metrics",
                        "#main/AnnotateReads/Valid_Reads",
                        "#main/AnnotateReads/Read1_error_rate",
                        "#main/AnnotateReads/Annotation_Read",
                        "#main/AnnotateReads/output",
                        "#main/AnnotateReads/validTcrReads",
                        "#main/AnnotateReads/validIgReads",
                        "#main/AnnotateReads/num_valid_tcr_reads",
                        "#main/AnnotateReads/num_valid_ig_reads"
                    ],
                    "requirements": [
                        {
                            "ramMin": 32000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AnnotateReads",
                    "in": [
                        {
                            "source": "#main/Internal_Settings/AbSeq_UMI",
                            "id": "#main/AnnotateReads/AbSeq_UMI"
                        },
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AnnotateReads/Extra_Seqs"
                        },
                        {
                            "source": "#main/QualityFilterOuter/Filter_Metrics",
                            "id": "#main/AnnotateReads/Filter_Metrics"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/AnnotateReads/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/AnnotateR1/Annotation_R1",
                            "id": "#main/AnnotateReads/R1_Annotation"
                        },
                        {
                            "source": "#main/AnnotateR1/R1_error_count_table",
                            "id": "#main/AnnotateReads/R1_error_count_table"
                        },
                        {
                            "source": "#main/AnnotateR1/R1_read_count_breakdown",
                            "id": "#main/AnnotateReads/R1_read_count_breakdown"
                        },
                        {
                            "source": "#main/AnnotateR2/Annot_R2",
                            "id": "#main/AnnotateReads/R2_Annotation"
                        },
                        {
                            "source": "#main/AnnotateR2/R2_Quality_Metrics",
                            "id": "#main/AnnotateReads/R2_Quality_Metrics"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/AnnotateReads/Run_Metadata"
                        },
                        {
                            "source": "#main/CheckReference/Target_Gene_Mapping",
                            "id": "#main/AnnotateReads/Target_Gene_Mapping"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/BundleLogs/logs_dir"
                    ],
                    "run": "#BundleLogs.cwl",
                    "id": "#main/BundleLogs",
                    "in": [
                        {
                            "source": [
                                "#main/AnnotateReads/output",
                                "#main/AnnotateR1/output",
                                "#main/AnnotateR2/output",
                                "#main/CheckReference/output",
                                "#main/GetDataTable/output",
                                "#main/Metrics/output",
                                "#main/AddtoBam/output",
                                "#main/AnnotateMolecules/output",
                                "#main/QualityFilterOuter/output",
                                "#main/CheckFastqs/log",
                                "#main/SplitAndSubsample/log",
                                "#main/MergeBAM/log",
                                "#main/Dense_to_Sparse_Datatable/output",
                                "#main/Dense_to_Sparse_Datatable_Unfiltered/output",
                                "#main/IndexBAM/log",
                                "#main/CellClassifier/log"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/BundleLogs/log_files"
                        }
                    ]
                },
                {
                    "run": "#Cell_Classifier.cwl",
                    "out": [
                        "#main/CellClassifier/cellTypePredictions",
                        "#main/CellClassifier/log"
                    ],
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/CellClassifier",
                    "in": [
                        {
                            "source": "#main/FindDataTableForCellClassifier/molsPerCellMatrixForCellClassifier",
                            "id": "#main/CellClassifier/molsPerCellMatrix"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/CheckFastqs/SubsampleSeed",
                        "#main/CheckFastqs/SubsamplingRatio",
                        "#main/CheckFastqs/FilesToSkipSplitAndSubsample",
                        "#main/CheckFastqs/FastqReadPairs",
                        "#main/CheckFastqs/Bead_Version",
                        "#main/CheckFastqs/Libraries",
                        "#main/CheckFastqs/ReadsList",
                        "#main/CheckFastqs/log"
                    ],
                    "run": "#CheckFastqs.cwl",
                    "id": "#main/CheckFastqs",
                    "in": [
                        {
                            "source": "#main/Internal_Settings/MinChunkSize",
                            "id": "#main/CheckFastqs/MinChunkSize"
                        },
                        {
                            "source": "#main/Reads",
                            "id": "#main/CheckFastqs/Reads"
                        },
                        {
                            "source": "#main/Subsample_Settings/Subsample_Reads",
                            "id": "#main/CheckFastqs/Subsample"
                        },
                        {
                            "source": "#main/Subsample_Settings/Subsample_Seed",
                            "id": "#main/CheckFastqs/Subsample_Seed"
                        }
                    ]
                },
                {
                    "run": "#CheckReference.cwl",
                    "out": [
                        "#main/CheckReference/Index",
                        "#main/CheckReference/Extra_Seqs",
                        "#main/CheckReference/Full_Genes",
                        "#main/CheckReference/output",
                        "#main/CheckReference/Transcript_Length",
                        "#main/CheckReference/GTF",
                        "#main/CheckReference/Target_Gene_Mapping"
                    ],
                    "requirements": [
                        {
                            "ramMin": 1000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/CheckReference",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "id": "#main/CheckReference/AbSeq_Reference"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/CheckReference/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/Reference",
                            "id": "#main/CheckReference/Reference"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/CheckReference/Run_Metadata"
                        }
                    ]
                },
                {
                    "run": "#DensetoSparse.cwl",
                    "scatter": [
                        "#main/Dense_to_Sparse_Datatable/Dense_Data_Table"
                    ],
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_File/Cell_Order",
                            "id": "#main/Dense_to_Sparse_Datatable/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Dense_Data_Tables",
                            "id": "#main/Dense_to_Sparse_Datatable/Dense_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Gene_List",
                            "id": "#main/Dense_to_Sparse_Datatable/Gene_List"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Dense_to_Sparse_Datatable/Run_Metadata"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 16000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/Dense_to_Sparse_Datatable",
                    "out": [
                        "#main/Dense_to_Sparse_Datatable/Data_Tables",
                        "#main/Dense_to_Sparse_Datatable/output"
                    ]
                },
                {
                    "run": "#DensetoSparse.cwl",
                    "scatter": [
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/Dense_Data_Table"
                    ],
                    "in": [
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Dense_Data_Tables_Unfiltered",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Dense_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Gene_List",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Gene_List"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Run_Metadata"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 16000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/Dense_to_Sparse_Datatable_Unfiltered",
                    "out": [
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/Data_Tables",
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/output"
                    ]
                },
                {
                    "out": [
                        "#main/Dense_to_Sparse_File/Cell_Order"
                    ],
                    "run": "#DensetoSparseFile.cwl",
                    "id": "#main/Dense_to_Sparse_File",
                    "in": [
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/Dense_to_Sparse_File/GDT_cell_order"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/FindDataTableForCellClassifier/molsPerCellMatrixForCellClassifier"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": {
                                    "items": "File",
                                    "type": "array"
                                },
                                "id": "#main/FindDataTableForCellClassifier/c174ddb5-9fdb-4dae-a1c5-b5666a631cc7/dataTables"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "type": "File",
                                "id": "#main/FindDataTableForCellClassifier/c174ddb5-9fdb-4dae-a1c5-b5666a631cc7/molsPerCellMatrixForCellClassifier"
                            }
                        ],
                        "id": "#main/FindDataTableForCellClassifier/c174ddb5-9fdb-4dae-a1c5-b5666a631cc7",
                        "expression": "${\n  for (var i = 0; i < inputs.dataTables.length; i++) {\n    var dataTable = inputs.dataTables[i];\n    if (dataTable.basename.indexOf(\"_RSEC_MolsPerCell.csv\") >= 0) {\n      return({molsPerCellMatrixForCellClassifier: dataTable});\n    }\n  }\n  return({molsPerCellMatrixForCellClassifier: null});\n}",
                        "class": "ExpressionTool"
                    },
                    "id": "#main/FindDataTableForCellClassifier",
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_Datatable/Data_Tables",
                            "id": "#main/FindDataTableForCellClassifier/dataTables"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/GetDataTable/Tag_Calls",
                        "#main/GetDataTable/Molecular_Annotation",
                        "#main/GetDataTable/Corrected_Molecular_Annotation",
                        "#main/GetDataTable/Tag_Annotation",
                        "#main/GetDataTable/Annot_Files",
                        "#main/GetDataTable/Cell_Label_Filter",
                        "#main/GetDataTable/Dense_Data_Tables",
                        "#main/GetDataTable/Dense_Data_Tables_Unfiltered",
                        "#main/GetDataTable/Expression_Data",
                        "#main/GetDataTable/Expression_Data_Unfiltered",
                        "#main/GetDataTable/Bioproduct_Stats",
                        "#main/GetDataTable/UMI_Adjusted_CellLabel_Stats",
                        "#main/GetDataTable/Putative_Cells_Origin",
                        "#main/GetDataTable/Protein_Aggregates_Experimental",
                        "#main/GetDataTable/Trueno_out",
                        "#main/GetDataTable/Trueno_zip",
                        "#main/GetDataTable/output",
                        "#main/GetDataTable/Cell_Order",
                        "#main/GetDataTable/Gene_List"
                    ],
                    "run": "#GetDataTable.cwl",
                    "id": "#main/GetDataTable",
                    "in": [
                        {
                            "source": "#main/CheckReference/Full_Genes",
                            "id": "#main/GetDataTable/Full_Genes"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Gene_Status_List",
                            "id": "#main/GetDataTable/Gene_Status_List"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Max_Count",
                            "id": "#main/GetDataTable/Max_Count"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Mol_Annot_List",
                            "id": "#main/GetDataTable/Molecule_Annotation_List"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/GetDataTable/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/GetDataTable/Run_Metadata"
                        },
                        {
                            "source": "#main/AnnotateReads/Seq_Metrics",
                            "id": "#main/GetDataTable/Seq_Metrics"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Tag_Sample_Names",
                            "id": "#main/GetDataTable/Tag_Names"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Total_Molecules",
                            "id": "#main/GetDataTable/Total_Molecules"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/IndexBAM/Index",
                        "#main/IndexBAM/log"
                    ],
                    "run": "#IndexBAM.cwl",
                    "id": "#main/IndexBAM",
                    "in": [
                        {
                            "source": "#main/MergeBAM/Final_Bam",
                            "id": "#main/IndexBAM/BamFile"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Internal_Settings/Read_Filter_Off",
                        "#main/Internal_Settings/Barcode_Num",
                        "#main/Internal_Settings/Seq_Run",
                        "#main/Internal_Settings/AbSeq_UMI",
                        "#main/Internal_Settings/Use_DBEC",
                        "#main/Internal_Settings/Extra_Seqs",
                        "#main/Internal_Settings/MinChunkSize",
                        "#main/Internal_Settings/NumRecordsPerSplit",
                        "#main/Internal_Settings/Target_analysis",
                        "#main/Internal_Settings/Subsample_Tags",
                        "#main/Internal_Settings/VDJ_VGene_Evalue",
                        "#main/Internal_Settings/VDJ_JGene_Evalue"
                    ],
                    "in": [],
                    "run": "#InternalSettings.cwl",
                    "id": "#main/Internal_Settings",
                    "label": "Internal Settings"
                },
                {
                    "out": [
                        "#main/MergeBAM/Final_Bam",
                        "#main/MergeBAM/log"
                    ],
                    "run": "#MergeBAM.cwl",
                    "id": "#main/MergeBAM",
                    "in": [
                        {
                            "source": "#main/AddtoBam/Annotated_Bam",
                            "id": "#main/MergeBAM/BamFiles"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Base_Name",
                            "id": "#main/MergeBAM/Run_Name"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/MergeBAM/Sample_Tags_Version"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/MergeMultiplex/Multiplex_out"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": {
                                    "items": [
                                        "null",
                                        "File"
                                    ],
                                    "type": "array"
                                },
                                "id": "#main/MergeMultiplex/8e7f752c-1505-4d65-81b3-f91fcd83b679/SampleTag_Files"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    {
                                        "items": "File",
                                        "type": "array"
                                    }
                                ],
                                "id": "#main/MergeMultiplex/8e7f752c-1505-4d65-81b3-f91fcd83b679/Multiplex_out"
                            }
                        ],
                        "id": "#main/MergeMultiplex/8e7f752c-1505-4d65-81b3-f91fcd83b679",
                        "expression": "${\n  var fp_array = [];\n  for (var i = 0; i < inputs.SampleTag_Files.length; i++) {\n    var fp = inputs.SampleTag_Files[i];\n    if (fp != null) {\n      fp_array.push(fp);\n    }\n  }\n  return({\"Multiplex_out\": fp_array});\n}",
                        "class": "ExpressionTool"
                    },
                    "id": "#main/MergeMultiplex",
                    "in": [
                        {
                            "source": [
                                "#main/GetDataTable/Trueno_out",
                                "#main/Metrics/Sample_Tag_Out"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/MergeMultiplex/SampleTag_Files"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Metadata_Settings/Run_Metadata",
                        "#main/Metadata_Settings/Run_Base_Name"
                    ],
                    "run": "#Metadata.cwl",
                    "id": "#main/Metadata_Settings",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "id": "#main/Metadata_Settings/AbSeq_Reference"
                        },
                        {
                            "valueFrom": "Targeted",
                            "id": "#main/Metadata_Settings/Assay"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Basic_Algo_Only",
                            "id": "#main/Metadata_Settings/Basic_Algo_Only"
                        },
                        {
                            "source": "#main/CheckFastqs/Bead_Version",
                            "id": "#main/Metadata_Settings/Bead_Version"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                            "id": "#main/Metadata_Settings/Exact_Cell_Count"
                        },
                        {
                            "source": "#main/CheckFastqs/Libraries",
                            "id": "#main/Metadata_Settings/Libraries"
                        },
                        {
                            "valueFrom": "BD Rhapsody Targeted Analysis Pipeline",
                            "id": "#main/Metadata_Settings/Pipeline_Name"
                        },
                        {
                            "source": "#main/Version/version",
                            "id": "#main/Metadata_Settings/Pipeline_Version"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                            "id": "#main/Metadata_Settings/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/CheckFastqs/ReadsList",
                            "id": "#main/Metadata_Settings/Reads"
                        },
                        {
                            "source": "#main/Reference",
                            "id": "#main/Metadata_Settings/Reference"
                        },
                        {
                            "source": "#main/Name_Settings/Run_Name",
                            "id": "#main/Metadata_Settings/Run_Name"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Tag_Sample_Names",
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
                            "source": "#main/Subsample_Settings/Subsample_Reads",
                            "id": "#main/Metadata_Settings/Subsample"
                        },
                        {
                            "source": "#main/Subsample_Settings/Subsample_Seed",
                            "id": "#main/Metadata_Settings/Subsample_Seed"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/Metadata_Settings/VDJ_Version"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Metrics/Metrics_Summary",
                        "#main/Metrics/Metrics_Archive",
                        "#main/Metrics/output",
                        "#main/Metrics/Sample_Tag_Out"
                    ],
                    "run": "#Metrics.cwl",
                    "id": "#main/Metrics",
                    "in": [
                        {
                            "source": "#main/GetDataTable/Annot_Files",
                            "id": "#main/Metrics/Annot_Files"
                        },
                        {
                            "source": "#main/AnnotateReads/Read1_error_rate",
                            "id": "#main/Metrics/Read1_error_rate"
                        },
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/Metrics/Run_Metadata"
                        },
                        {
                            "source": "#main/GetDataTable/Trueno_zip",
                            "id": "#main/Metrics/Sample_Tag_Archives"
                        },
                        {
                            "source": "#main/Internal_Settings/Seq_Run",
                            "id": "#main/Metrics/Seq_Run"
                        },
                        {
                            "source": "#main/GetDataTable/UMI_Adjusted_CellLabel_Stats",
                            "id": "#main/Metrics/UMI_Adjusted_Stats"
                        },
                        {
                            "source": "#main/VDJ_Compile_Results/vdjMetricsJson",
                            "id": "#main/Metrics/vdjMetricsJson"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Multiplexing_Settings/Tag_Sample_Names",
                        "#main/Multiplexing_Settings/Sample_Tags_Version"
                    ],
                    "in": [
                        {
                            "source": "#main/Sample_Tags_Version",
                            "id": "#main/Multiplexing_Settings/_Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Tag_Names",
                            "id": "#main/Multiplexing_Settings/_Tag_Sample_Names"
                        }
                    ],
                    "run": "#MultiplexingSettings.cwl",
                    "id": "#main/Multiplexing_Settings",
                    "label": "Multiplexing Settings"
                },
                {
                    "out": [
                        "#main/Name_Settings/Run_Name"
                    ],
                    "in": [
                        {
                            "source": "#main/Run_Name",
                            "id": "#main/Name_Settings/_Run_Name"
                        }
                    ],
                    "run": "#NameSettings.cwl",
                    "id": "#main/Name_Settings",
                    "label": "Name Settings"
                },
                {
                    "out": [
                        "#main/PairReadFiles/ReadPairs"
                    ],
                    "run": "#PairReadFiles.cwl",
                    "id": "#main/PairReadFiles",
                    "in": [
                        {
                            "source": "#main/CheckFastqs/FastqReadPairs",
                            "id": "#main/PairReadFiles/FastqReadPairs"
                        },
                        {
                            "source": "#main/SplitAndSubsample/SplitAndSubsampledFastqs",
                            "id": "#main/PairReadFiles/Reads"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Putative_Cell_Calling_Settings/Putative_Cell_Call",
                        "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                        "#main/Putative_Cell_Calling_Settings/Basic_Algo_Only"
                    ],
                    "in": [
                        {
                            "source": "#main/Basic_Algo_Only",
                            "id": "#main/Putative_Cell_Calling_Settings/_Basic_Algo_Only"
                        },
                        {
                            "source": "#main/Exact_Cell_Count",
                            "id": "#main/Putative_Cell_Calling_Settings/_Exact_Cell_Count"
                        },
                        {
                            "source": "#main/Putative_Cell_Call",
                            "id": "#main/Putative_Cell_Calling_Settings/_Putative_Cell_Call"
                        }
                    ],
                    "run": "#PutativeCellSettings.cwl",
                    "id": "#main/Putative_Cell_Calling_Settings",
                    "label": "Putative Cell Calling Settings"
                },
                {
                    "out": [
                        "#main/QualityFilterOuter/Filter_Metrics",
                        "#main/QualityFilterOuter/R1",
                        "#main/QualityFilterOuter/R2",
                        "#main/QualityFilterOuter/output"
                    ],
                    "run": "#QualityFilterOuter.cwl",
                    "id": "#main/QualityFilterOuter",
                    "in": [
                        {
                            "source": "#main/Metadata_Settings/Run_Metadata",
                            "id": "#main/QualityFilterOuter/Run_Metadata"
                        },
                        {
                            "source": "#main/PairReadFiles/ReadPairs",
                            "id": "#main/QualityFilterOuter/Split_Read_Pairs"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/SplitAndSubsample/SplitAndSubsampledFastqs",
                        "#main/SplitAndSubsample/log"
                    ],
                    "run": "#SplitAndSubsample.cwl",
                    "id": "#main/SplitAndSubsample",
                    "in": [
                        {
                            "source": "#main/Reads",
                            "id": "#main/SplitAndSubsample/Fastqs"
                        },
                        {
                            "source": "#main/CheckFastqs/FilesToSkipSplitAndSubsample",
                            "id": "#main/SplitAndSubsample/FilesToSkipSplitAndSubsample"
                        },
                        {
                            "source": "#main/Internal_Settings/NumRecordsPerSplit",
                            "id": "#main/SplitAndSubsample/NumRecordsPerSplit"
                        },
                        {
                            "source": "#main/CheckFastqs/SubsamplingRatio",
                            "id": "#main/SplitAndSubsample/SubsampleRatio"
                        },
                        {
                            "source": "#main/CheckFastqs/SubsampleSeed",
                            "id": "#main/SplitAndSubsample/SubsampleSeed"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Start_Time/Start_Time"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "type": "string",
                                "id": "#main/Start_Time/dc4e9fd7-92dc-4aca-80ad-76601aaaf6ad/Start_Time"
                            }
                        ],
                        "id": "#main/Start_Time/dc4e9fd7-92dc-4aca-80ad-76601aaaf6ad",
                        "expression": "${   \n  var today = new Date();\n  var date = today.toString()\n  return ({Start_Time: date});\n} ",
                        "class": "ExpressionTool"
                    },
                    "id": "#main/Start_Time",
                    "in": []
                },
                {
                    "out": [
                        "#main/Subsample_Settings/Subsample_Reads",
                        "#main/Subsample_Settings/Subsample_Seed"
                    ],
                    "in": [
                        {
                            "source": "#main/Subsample",
                            "id": "#main/Subsample_Settings/_Subsample_Reads"
                        },
                        {
                            "source": "#main/Subsample_seed",
                            "id": "#main/Subsample_Settings/_Subsample_Seed"
                        }
                    ],
                    "run": "#SubsampleSettings.cwl",
                    "id": "#main/Subsample_Settings",
                    "label": "Subsample Settings"
                },
                {
                    "out": [
                        "#main/Uncompress_Datatables/Uncompressed_Data_Tables",
                        "#main/Uncompress_Datatables/Uncompressed_Expression_Matrix"
                    ],
                    "run": "#UncompressDatatables.cwl",
                    "id": "#main/Uncompress_Datatables",
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_Datatable/Data_Tables",
                            "id": "#main/Uncompress_Datatables/Compressed_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Expression_Data",
                            "id": "#main/Uncompress_Datatables/Compressed_Expression_Matrix"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_Assemble_and_Annotate_Contigs_IG/igCalls"
                    ],
                    "run": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl",
                    "id": "#main/VDJ_Assemble_and_Annotate_Contigs_IG",
                    "in": [
                        {
                            "source": "#main/VDJ_Preprocess_Reads_IG/RSEC_Reads_Fastq",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_IG/RSEC_Reads_Fastq"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_IG/VDJ_Version"
                        },
                        {
                            "source": "#main/VDJ_Preprocess_Reads_IG/num_cores",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_IG/num_cores"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/tcrCalls"
                    ],
                    "run": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl",
                    "id": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR",
                    "in": [
                        {
                            "source": "#main/VDJ_Preprocess_Reads_TCR/RSEC_Reads_Fastq",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/RSEC_Reads_Fastq"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/VDJ_Version"
                        },
                        {
                            "source": "#main/VDJ_Preprocess_Reads_TCR/num_cores",
                            "id": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/num_cores"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_Compile_Results/vdjCellsDatatable",
                        "#main/VDJ_Compile_Results/vdjCellsDatatableUncorrected",
                        "#main/VDJ_Compile_Results/vdjDominantContigs",
                        "#main/VDJ_Compile_Results/vdjUnfilteredContigs",
                        "#main/VDJ_Compile_Results/vdjMetricsJson",
                        "#main/VDJ_Compile_Results/vdjMetricsCsv",
                        "#main/VDJ_Compile_Results/vdjReadsPerCellByChainTypeFigure"
                    ],
                    "run": "#VDJ_Compile_Results.cwl",
                    "id": "#main/VDJ_Compile_Results",
                    "in": [
                        {
                            "source": "#main/AnnotateReads/Seq_Metrics",
                            "id": "#main/VDJ_Compile_Results/Seq_Metrics"
                        },
                        {
                            "source": "#main/CellClassifier/cellTypePredictions",
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
                            "source": "#main/VDJ_GatherIGCalls/gatheredCalls",
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
                            "source": "#main/VDJ_GatherTCRCalls/gatheredCalls",
                            "id": "#main/VDJ_Compile_Results/tcrCalls"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_Compile_Results/vdjVersion"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_GatherIGCalls/gatheredCalls"
                    ],
                    "run": "#VDJ_GatherCalls.cwl",
                    "id": "#main/VDJ_GatherIGCalls",
                    "in": [
                        {
                            "source": "#main/VDJ_Assemble_and_Annotate_Contigs_IG/igCalls",
                            "id": "#main/VDJ_GatherIGCalls/theCalls"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_GatherTCRCalls/gatheredCalls"
                    ],
                    "run": "#VDJ_GatherCalls.cwl",
                    "id": "#main/VDJ_GatherTCRCalls",
                    "in": [
                        {
                            "source": "#main/VDJ_Assemble_and_Annotate_Contigs_TCR/tcrCalls",
                            "id": "#main/VDJ_GatherTCRCalls/theCalls"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_Preprocess_Reads_IG/RSEC_Reads_Fastq",
                        "#main/VDJ_Preprocess_Reads_IG/num_splits",
                        "#main/VDJ_Preprocess_Reads_IG/num_cores"
                    ],
                    "run": "#VDJ_Preprocess_Reads.cwl",
                    "id": "#main/VDJ_Preprocess_Reads_IG",
                    "in": [
                        {
                            "source": "#main/AnnotateReads/validIgReads",
                            "id": "#main/VDJ_Preprocess_Reads_IG/Valid_Reads_Fastq"
                        },
                        {
                            "source": "#main/AnnotateReads/num_valid_ig_reads",
                            "id": "#main/VDJ_Preprocess_Reads_IG/num_valid_reads"
                        },
                        {
                            "valueFrom": "BCR",
                            "id": "#main/VDJ_Preprocess_Reads_IG/vdj_type"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_Preprocess_Reads_TCR/RSEC_Reads_Fastq",
                        "#main/VDJ_Preprocess_Reads_TCR/num_splits",
                        "#main/VDJ_Preprocess_Reads_TCR/num_cores"
                    ],
                    "run": "#VDJ_Preprocess_Reads.cwl",
                    "id": "#main/VDJ_Preprocess_Reads_TCR",
                    "in": [
                        {
                            "source": "#main/AnnotateReads/validTcrReads",
                            "id": "#main/VDJ_Preprocess_Reads_TCR/Valid_Reads_Fastq"
                        },
                        {
                            "source": "#main/AnnotateReads/num_valid_tcr_reads",
                            "id": "#main/VDJ_Preprocess_Reads_TCR/num_valid_reads"
                        },
                        {
                            "valueFrom": "TCR",
                            "id": "#main/VDJ_Preprocess_Reads_TCR/vdj_type"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_Settings/VDJ_Version"
                    ],
                    "in": [
                        {
                            "source": "#main/VDJ_Version",
                            "id": "#main/VDJ_Settings/_VDJ_Version"
                        }
                    ],
                    "run": "#VDJ_Settings.cwl",
                    "id": "#main/VDJ_Settings",
                    "label": "VDJ Settings"
                },
                {
                    "out": [
                        "#main/Version/version"
                    ],
                    "run": "#Version.cwl",
                    "id": "#main/Version",
                    "in": []
                }
            ],
            "outputs": [
                {
                    "outputSource": "#main/GetDataTable/Bioproduct_Stats",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Bioproduct_Stats",
                    "label": "Bioproduct Statistics"
                },
                {
                    "outputSource": "#main/GetDataTable/Cell_Label_Filter",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Cell_Label_Filter",
                    "label": "Cell Label Filter"
                },
                {
                    "outputSource": "#main/Uncompress_Datatables/Uncompressed_Data_Tables",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Data_Tables",
                    "label": "Data Tables"
                },
                {
                    "outputSource": "#main/Dense_to_Sparse_Datatable_Unfiltered/Data_Tables",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Data_Tables_Unfiltered",
                    "label": "Unfiltered Data Tables"
                },
                {
                    "outputSource": "#main/Uncompress_Datatables/Uncompressed_Expression_Matrix",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Expression_Data",
                    "label": "Expression Matrix"
                },
                {
                    "outputSource": "#main/GetDataTable/Expression_Data_Unfiltered",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Expression_Data_Unfiltered",
                    "label": "Unfiltered Expression Matrix"
                },
                {
                    "outputSource": "#main/MergeBAM/Final_Bam",
                    "type": "File",
                    "id": "#main/Final_Bam",
                    "label": "Final BAM File"
                },
                {
                    "outputSource": "#main/IndexBAM/Index",
                    "type": "File",
                    "id": "#main/Final_Bam_Index",
                    "label": "Final BAM Index"
                },
                {
                    "outputSource": "#main/CellClassifier/cellTypePredictions",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/ImmuneCellClassification(Experimental)",
                    "label": "Immune Cell Classification (Experimental)"
                },
                {
                    "outputSource": "#main/BundleLogs/logs_dir",
                    "type": "Directory",
                    "id": "#main/Logs",
                    "label": "Pipeline Logs"
                },
                {
                    "outputSource": "#main/Metrics/Metrics_Summary",
                    "type": "File",
                    "id": "#main/Metrics_Summary",
                    "label": "Metrics Summary"
                },
                {
                    "outputSource": "#main/MergeMultiplex/Multiplex_out",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Multiplex"
                },
                {
                    "outputSource": "#main/GetDataTable/Protein_Aggregates_Experimental",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Protein_Aggregates_Experimental",
                    "label": "Protein Aggregates (Experimental)"
                },
                {
                    "outputSource": "#main/GetDataTable/Putative_Cells_Origin",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Putative_Cells_Origin",
                    "label": "Putative Cells Origin"
                },
                {
                    "outputSource": "#main/VDJ_Compile_Results/vdjCellsDatatable",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/vdjCellsDatatable",
                    "label": "vdjCellsDatatable"
                },
                {
                    "outputSource": "#main/VDJ_Compile_Results/vdjCellsDatatableUncorrected",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/vdjCellsDatatableUncorrected",
                    "label": "vdjCellsDatatableUncorrected"
                },
                {
                    "outputSource": "#main/VDJ_Compile_Results/vdjDominantContigs",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/vdjDominantContigs",
                    "label": "vdjDominantContigs"
                },
                {
                    "outputSource": "#main/VDJ_Compile_Results/vdjMetricsCsv",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/vdjMetricsCsv",
                    "label": "vdjMetricsCsv"
                },
                {
                    "outputSource": "#main/VDJ_Compile_Results/vdjUnfilteredContigs",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/vdjUnfilteredContigs",
                    "label": "vdjUnfilteredContigs"
                }
            ],
            "id": "#main",
            "class": "Workflow"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#MergeBAM.cwl/BamFiles"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MergeBAM.cwl/Run_Name"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MergeBAM.cwl/Sample_Tags_Version"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "samtools_merge.log",
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_final.BAM"
                    },
                    "type": "File",
                    "id": "#MergeBAM.cwl/Final_Bam"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#MergeBAM.cwl/log"
                }
            ],
            "baseCommand": [
                "samtools",
                "merge"
            ],
            "id": "#MergeBAM.cwl",
            "arguments": [
                {
                    "prefix": "-@",
                    "valueFrom": "$(runtime.cores)"
                },
                {
                    "position": 0,
                    "valueFrom": "${\n    if (inputs.Sample_Tags_Version) {\n        return \"Combined_\" + inputs.Run_Name + \"_final.BAM\"\n    } else {\n        return inputs.Run_Name + \"_final.BAM\"\n    }\n}"
                }
            ],
            "class": "CommandLineTool",
            "hints": [
                {
                    "coresMin": 4,
                    "class": "ResourceRequirement"
                }
            ]
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#Metadata.cwl/AbSeq_Reference"
                },
                {
                    "type": "string",
                    "id": "#Metadata.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#Metadata.cwl/Basic_Algo_Only"
                },
                {
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "string",
                                    "name": "#Metadata.cwl/Bead_Version/Library"
                                },
                                {
                                    "type": "string",
                                    "name": "#Metadata.cwl/Bead_Version/bead_version"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#Metadata.cwl/Bead_Version"
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
                        "int"
                    ],
                    "id": "#Metadata.cwl/Label_Version"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Libraries"
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
                        "int"
                    ],
                    "id": "#Metadata.cwl/Putative_Cell_Call"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#Metadata.cwl/Read_Filter_Off"
                },
                {
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#Metadata.cwl/Reads"
                },
                {
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#Metadata.cwl/Reference"
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
                            "items": "string",
                            "type": "array"
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
                        "float"
                    ],
                    "id": "#Metadata.cwl/Subsample"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#Metadata.cwl/Subsample_Seed"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#Metadata.cwl/Subsample_Tags"
                },
                {
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#Metadata.cwl/Supplemental_Reference"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/VDJ_Version"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "run_metadata.json",
            "outputs": [
                {
                    "outputBinding": {
                        "outputEval": "${  \n  var name = inputs.Run_Name;\n  if (name == null){\n    var libraries = inputs.Libraries;\n    name = libraries.split(',')[0];\n  }   \n  return(name)\n}   \n"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metadata.cwl/Run_Base_Name"
                },
                {
                    "type": "stdout",
                    "id": "#Metadata.cwl/Run_Metadata"
                }
            ],
            "baseCommand": "echo",
            "id": "#Metadata.cwl",
            "arguments": [
                {
                    "prefix": ""
                },
                {
                    "shellQuote": true,
                    "valueFrom": "${\n  var metadata = inputs;\n  var all_bv = {};\n  var customer_bv = \"Original (V1)\";\n  for (var i = 0; i < inputs.Bead_Version.length; i++) {\n      var BeadVer = inputs.Bead_Version[i];\n      var Library = BeadVer[\"Library\"];\n      var bead_version = BeadVer[\"bead_version\"];\n      all_bv[Library] = bead_version  \n      var short_bv =  bead_version.substring(0, 2);\n      if (short_bv == \"V2\"){\n        var customer_bv = \"Enhanced (V2)\";\n      }\n  }\n  metadata[\"Bead_Version\"] = all_bv;\n\n  var pipeline_name = inputs.Pipeline_Name;\n  var assay = inputs.Assay;\n  var version = inputs.Pipeline_Version;\n  var time = inputs.Start_Time;\n  var libraries = inputs.Libraries.split(\",\");\n  var i = 0;\n  var reference_list = []\n  if(inputs.Reference != null){\n      reference_list = reference_list.concat(inputs.Reference);\n  }\n  if(inputs.AbSeq_Reference != null){\n      reference_list = reference_list.concat(inputs.AbSeq_Reference);\n  }\n\n  var supplemental = \"\"\n  if(inputs.Supplemental_Reference != null){\n      supplemental = \"; Supplemental_Reference - \" + inputs.Supplemental_Reference[0][\"basename\"];\n  }\n  var references = [];\n  for (i = 0; i< reference_list.length; i++) {\n      if(reference_list[i] != null){\n          references.push(reference_list[i][\"basename\"]);\n      }\n  }\n  var parameters = [];\n  if(inputs.Sample_Tags_Version != null){\n      var tags = \"Sample Tag Version: \" + inputs.Sample_Tags_Version;\n  } else{ \n      var tags = \"Sample Tag Version: None\";\n  }\n  parameters.push(tags);\n\n  if(inputs.Sample_Tag_Names != null){\n      var tag_names = inputs.Sample_Tag_Names.join(\" ; \")\n      var tag_list = \"Sample Tag Names: \" + tag_names;\n  } else{\n      var tag_list = \"Sample Tag Names: None\";\n  }\n  parameters.push(tag_list);\n \n  if(inputs.VDJ_Version != null){\n      var vdj = \"VDJ Version: \" + inputs.VDJ_Version;\n  } else{ \n      var vdj = \"VDJ Version: None\";\n  }\n  parameters.push(vdj)\n\n  if(inputs.Subsample != null){\n      var subsample = \"Subsample: \" + inputs.Subsample;\n  } else{ \n      var subsample = \"Subsample: None\";\n  }   \n  parameters.push(subsample);\n\n  if(inputs.Putative_Cell_Call == 1){\n      var call = \"Putative Cell Calling Type: AbSeq\";\n  } else{ \n      var call = \"Putative Cell Calling Type: mRNA\";\n  }   \n  parameters.push(call)\n\n  if(inputs.Basic_Algo_Only){\n      var basic = \"Refined Putative Cell Calling: Off\";\n  } else{ \n      var basic = \"Refined Putative Cell Calling: On\";\n  }   \n  parameters.push(basic)\n\n  if(inputs.Exact_Cell_Count != null){\n      var cells = \"Exact Cell Count: \" + inputs.Exact_Cell_Count;\n  } else{ \n      var cells = \"Exact Cell Count: None\";\n  }   \n  parameters.push(cells)\n\n  var name = inputs.Run_Name;\n  if (name == null){\n    var libraries = inputs.Libraries.split(',');\n    name = libraries[0];\n  }        \n\n  var header = [\"####################\"];\n  header.push(\"## \" + pipeline_name + \" Version \" + version);\n  header.push(\"## Analysis Date - \" + time);\n  header.push(\"## Libraries - \" + libraries.join(' | ') + \" - Bead version detected: \" + customer_bv);\n  header.push(\"## References - \" + references.join(' | ') + supplemental);\n  header.push(\"## Parameters - \" + parameters.join(' | '));\n  header.push(\"####################\");\n  metadata[\"Output_Header\"] = header;\n  metadata[\"Run_Base_Name\"] = name;\n  var metadata_json = JSON.stringify(metadata);\n  return metadata_json;\n}\n"
                }
            ],
            "class": "CommandLineTool"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--annot-files"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Annot_Files"
                },
                {
                    "inputBinding": {
                        "prefix": "--read1-error-rate"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Read1_error_rate"
                },
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Run_Metadata"
                },
                {
                    "inputBinding": {
                        "prefix": "--sample-tag-archives",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#Metrics.cwl/Sample_Tag_Archives"
                },
                {
                    "inputBinding": {
                        "prefix": "--seq-run"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metrics.cwl/Seq_Run"
                },
                {
                    "inputBinding": {
                        "prefix": "--umi-adjusted-stats"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Metrics.cwl/UMI_Adjusted_Stats"
                },
                {
                    "inputBinding": {
                        "prefix": "--vdj-metrics-fp"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Metrics.cwl/vdjMetricsJson"
                }
            ],
            "requirements": [
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "internal-metrics-archive.tar.gz"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Metrics_Archive"
                },
                {
                    "outputBinding": {
                        "glob": "*_Metrics_Summary.csv"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Metrics_Summary"
                },
                {
                    "outputBinding": {
                        "glob": "*.zip"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#Metrics.cwl/Sample_Tag_Out"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_metrics.py"
            ],
            "class": "CommandLineTool",
            "id": "#Metrics.cwl"
        },
        {
            "inputs": [
                {
                    "default": "Targeted",
                    "type": "string",
                    "id": "#MultiplexingSettings.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#MultiplexingSettings.cwl/_Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/_Tag_Sample_Names"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MultiplexingSettings.cwl/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/Tag_Sample_Names"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  var enumifiedSampleTagsVersion = null;\n  if (inputs._Sample_Tags_Version) {\n  var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();\n  if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')\n  {\n    enumifiedSampleTagsVersion = 'hs';\n  }\n  else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')\n  {\n    enumifiedSampleTagsVersion = 'mm';\n  }\n  else if (_Sample_Tags_Version === 'no multiplexing')\n  {\n    enumifiedSampleTagsVersion = null;\n  }\n  else\n  {\n    throw new Error(\"Cannot parse Sample Tag Version: \" + inputs._Sample_Tags_Version);\n  }\n  }\n  var listTagNames = inputs._Tag_Sample_Names\n  var newTagNames = []\n  for (var num in listTagNames) {\n    var tag = listTagNames[num].replace(/[^A-Za-z0-9-+]/g,\"_\");\n    newTagNames.push(tag);    \n  }  \n  return ({\n  Tag_Sample_Names: newTagNames,\n  Sample_Tags_Version: enumifiedSampleTagsVersion\n  });\n}",
            "id": "#MultiplexingSettings.cwl"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#NameSettings.cwl/_Run_Name"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#NameSettings.cwl/Run_Name"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${ var name = inputs._Run_Name;\n   if (name != null) {\n     name = name.replace(/[\\W_]+/g,\"-\");}\n   return({'Run_Name' : name });\n }  ",
            "id": "#NameSettings.cwl"
        },
        {
            "inputs": [
                {
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/filename"
                                },
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/readFlag"
                                },
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/readPairId"
                                },
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/library"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#PairReadFiles.cwl/FastqReadPairs"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#PairReadFiles.cwl/Reads"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "PairReadFiles takes an array of split files and pairs them, such that an R1 file is transferred to the QualityFilter with its corresponding R2 file.\nThe original FASTQ files are paired in CheckFastqs and then split and sub-sampled in SplitAndSubsample. The pairing information is taken from CheckFastqs.\n",
            "id": "#PairReadFiles.cwl",
            "outputs": [
                {
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "File",
                                    "name": "#PairReadFiles.cwl/ReadPairs/R1"
                                },
                                {
                                    "type": "File",
                                    "name": "#PairReadFiles.cwl/ReadPairs/R2"
                                },
                                {
                                    "type": "int",
                                    "name": "#PairReadFiles.cwl/ReadPairs/readPairId"
                                },
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/ReadPairs/library"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#PairReadFiles.cwl/ReadPairs"
                }
            ],
            "expression": "${\n  // use the CheckFastqs read pairing information to create a dictionary\n  // using the original fastq file name without the extension as the key\n  var fastqReadPairs = {}\n  for (var i = 0; i < inputs.FastqReadPairs.length; i++) {\n    var fileDict = inputs.FastqReadPairs[i];\n    var filename = fileDict[\"filename\"];\n\n    if (!fastqReadPairs[filename]) {\n      fastqReadPairs[filename] = {\n        readPairId: null,\n        readFlag: null,\n        library: null,\n      };\n    }\n    else {\n      throw new Error(\"Found non-unique fastq filename '\" + filename + \"' in the FastqReadPairs dictionary from CheckFastqs.\")\n    }\n\n    fastqReadPairs[filename].readPairId = fileDict[\"readPairId\"]\n    fastqReadPairs[filename].readFlag = fileDict[\"readFlag\"]\n    fastqReadPairs[filename].library = fileDict[\"library\"]\n  }\n\n  // now loop through the input read files which could\n  // be the original fastq files if no sub-sampling has\n  // been done, or the sub-sampled fastq files\n  var readPairs = {}\n  for (var i = 0; i < inputs.Reads.length; i++) {\n\n    // Set the fileDict to null\n    var fileDict = null;\n\n    // Get the fastq file\n    var fastqFile = inputs.Reads[i];\n\n    // Remove the .gz from the end of the filename\n    var fileNoGzExt = fastqFile.basename.replace(/.gz$/i, \"\");\n\n    // Remove the next file extension if it exists\n    var fileArrayWithExt = fileNoGzExt.split(\".\");\n    // If an extension exists, splice the array\n    var fileArrayNoExt = null;\n    if (fileArrayWithExt.length > 1) {\n      fileArrayNoExt = fileArrayWithExt.splice(0, fileArrayWithExt.length-1);\n    } else {\n      // No file extension exists, so use the whole array\n      fileArrayNoExt = fileArrayWithExt\n    }\n    var fileRootname = fileArrayNoExt.join(\".\")\n\n    // if the original files were sub-sampled\n    // get the original file and the chunk id\n    if (fileRootname.indexOf(\"-\") != -1) {\n      // Split on the dash to get the name of\n      // the original file and the chunk id\n      // The original file name can also have dashes\n      var chunkFileArray = fileRootname.split(\"-\");\n\n      // Get the original file rootname and chunk id\n      // The rootname without the chunk id and file\n      // extension is the key from CheckFastqs\n      // The chunk id is used later to create a new unique\n      // read pair id for all sub-sampled fastq files\n\n      // The rootname array should contain all elements up to the last dash\n      var fileRootnameArray = chunkFileArray.splice(0, chunkFileArray.length-1);\n      var fileRootnameNoChunkId = fileRootnameArray.join(\"-\");\n\n      // The chunk id is the last element in the array\n      // representing the content after the last dash\n      var orgChunkId = chunkFileArray.pop();\n\n      // if there is no chunk id, use an arbitrary number\n      // the chunk id is unique when the files are sub-sampled\n      // and does not need to be unique when the files are not sub-sampled\n      var chunkId = 9999;\n      if (orgChunkId) {\n        // cast to an integer\n        chunkId = parseInt(orgChunkId);\n      }\n      // double check that we have a chunk id\n      if (chunkId === undefined || chunkId === null) {\n        throw new Error(\"The fastq file sub-sampling id could not be determined!\");\n      }\n\n      // The file rootname without the chunk id and file extension\n      // should match the original file rootname from CheckFastqs\n      // The original file rootname from CheckFastqs is the key for\n      // the dictionary containing the original unique pair id\n      var fileDict = fastqReadPairs[fileRootnameNoChunkId];\n    }\n\n    // If the files are not sub-sampled or the fileDict\n    // is not found, then try to use the original\n    // file rootname without the file extension as the key\n    if (fileDict === undefined || fileDict === null) {\n\n      // if the original files were not sub-sampled,\n      // use the original file rootname and an arbitrary chunk id\n      var chunkId = 9999;\n\n      var fileDict = fastqReadPairs[fileRootname];\n\n      // If the fileDict for this file rootname is not found,\n      // then the filenames are in an unexpected format and\n      // the code to parse the filenames in CheckFastqs,\n      // SplitAndSubsample and here need to match\n      if (fileDict === undefined || fileDict === null) {\n        // Create an error\n        if (fileDict === undefined || fileDict === null) {\n          throw new Error(\"Cannot find the fastq read pair information for '\" + fastqFile.basename + \"'.\");\n        }\n      }\n    }\n\n    // Get the pairing information from CheckFastqs\n    var readPairId = fileDict[\"readPairId\"];\n    var library = fileDict[\"library\"];\n    var flag = fileDict[\"readFlag\"];\n\n    // Add the chunkId to create a new unique read pair id\n    // for each file (sub-sampled or not)\n    var chunkReadPairId = readPairId + \"_\" + chunkId;\n\n    // Create a dictionary for each pair of files\n    if (!readPairs[chunkReadPairId]) {\n      readPairs[chunkReadPairId] = {\n        R1: null,\n        R2: null,\n        library: library,\n        readPairId: null,\n      };\n    }\n    // add in the R1 and R2 files, depending on the flag\n    if (flag === \"R1\") {\n      readPairs[chunkReadPairId].R1 = fastqFile\n    } else if (flag === \"R2\") {\n      readPairs[chunkReadPairId].R2 = fastqFile\n    }\n  }\n  // we are not interested in the read pair ids in readPairs\n  // flatten into an array of objects\n  var readPairsList = [];\n  var i = 1;\n  for (var key in readPairs) {\n    if (readPairs.hasOwnProperty(key)) {\n      var readPair = readPairs[key];\n      readPair.readPairId = i;\n      readPairsList.push(readPair);\n      i++;\n    }\n  }\n  // pass this array to the record array named \"ReadPairs\" on the CWL layer\n  return {ReadPairs: readPairsList}\n}",
            "class": "ExpressionTool"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Basic_Algo_Only"
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
                        "Any"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Putative_Cell_Call"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#PutativeCellSettings.cwl/Basic_Algo_Only"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/Exact_Cell_Count"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/Putative_Cell_Call"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  // the basic algorithm flag defaults to false\n  var basicAlgOnlyFlag = false;\n  // the user can set the basic algorithm flag\n  if (inputs._Basic_Algo_Only) {\n    basicAlgOnlyFlag = inputs._Basic_Algo_Only;\n  }\n  // convert the Putative_Cell_Call from a string to an integer\n  var putativeCellCallInt = 0;\n  if (inputs._Putative_Cell_Call) {\n    if (inputs._Putative_Cell_Call === \"mRNA\") {\n      putativeCellCallInt = 0;\n    }\n    else if (inputs._Putative_Cell_Call == \"AbSeq_Experimental\" || inputs._Putative_Cell_Call == \"AbSeq (Experimental)\") {\n      putativeCellCallInt = 1;\n      // for protein-only cell calling, we only have the basic algorithm\n      basicAlgOnlyFlag = true;\n    }\n    else if (inputs._Putative_Cell_Call == \"mRNA_and_AbSeq\") {\n      putativeCellCallInt = 2;\n    }\n  }\n  // check the exact cell count\n  if (inputs._Exact_Cell_Count) {\n    if (inputs._Exact_Cell_Count < 1) {\n      throw(\"Illogical value for exact cell count: \" + inputs._Exact_Cell_Count);\n    }\n  }\n  return ({\n    Putative_Cell_Call: putativeCellCallInt,\n    Exact_Cell_Count: inputs._Exact_Cell_Count,\n    Basic_Algo_Only: basicAlgOnlyFlag,\n  });\n}",
            "id": "#PutativeCellSettings.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--run-metadata"
                    },
                    "type": "File",
                    "id": "#QualityFilter.cwl/Run_Metadata"
                },
                {
                    "type": {
                        "fields": [
                            {
                                "inputBinding": {
                                    "prefix": "--r1"
                                },
                                "type": "File",
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/R1"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--r2"
                                },
                                "type": "File",
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/R2"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--read-pair-id"
                                },
                                "type": "int",
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/readPairId"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--library"
                                },
                                "type": "string",
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/library"
                            }
                        ],
                        "type": "record"
                    },
                    "id": "#QualityFilter.cwl/Split_Read_Pairs"
                }
            ],
            "requirements": [
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*read_quality.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualityFilter.cwl/Filter_Metrics"
                },
                {
                    "outputBinding": {
                        "glob": "*_R1*.fastq.gz"
                    },
                    "type": "File",
                    "id": "#QualityFilter.cwl/R1"
                },
                {
                    "outputBinding": {
                        "glob": "*_R2*.fastq.gz"
                    },
                    "type": "File",
                    "id": "#QualityFilter.cwl/R2"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#QualityFilter.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_quality_filter.py"
            ],
            "class": "CommandLineTool",
            "id": "#QualityFilter.cwl"
        },
        {
            "inputs": [
                {
                    "type": "File",
                    "id": "#QualityFilterOuter.cwl/Run_Metadata"
                },
                {
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "File",
                                    "name": "#QualityFilterOuter.cwl/Split_Read_Pairs/R1"
                                },
                                {
                                    "type": "File",
                                    "name": "#QualityFilterOuter.cwl/Split_Read_Pairs/R2"
                                },
                                {
                                    "type": "int",
                                    "name": "#QualityFilterOuter.cwl/Split_Read_Pairs/readPairId"
                                },
                                {
                                    "type": "string",
                                    "name": "#QualityFilterOuter.cwl/Split_Read_Pairs/library"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#QualityFilterOuter.cwl/Split_Read_Pairs"
                }
            ],
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
            "outputs": [
                {
                    "outputSource": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Filter_Metrics",
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#QualityFilterOuter.cwl/Filter_Metrics"
                },
                {
                    "outputSource": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/R1",
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#QualityFilterOuter.cwl/R1"
                },
                {
                    "outputSource": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/R2",
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#QualityFilterOuter.cwl/R2"
                },
                {
                    "outputSource": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/output",
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#QualityFilterOuter.cwl/output"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "scatter": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Split_Read_Pairs",
                    "out": [
                        "#QualityFilterOuter.cwl/Quality_Filter_Scatter/R1",
                        "#QualityFilterOuter.cwl/Quality_Filter_Scatter/R2",
                        "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Filter_Metrics",
                        "#QualityFilterOuter.cwl/Quality_Filter_Scatter/output"
                    ],
                    "run": "#QualityFilter.cwl",
                    "id": "#QualityFilterOuter.cwl/Quality_Filter_Scatter",
                    "in": [
                        {
                            "source": "#QualityFilterOuter.cwl/Run_Metadata",
                            "id": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Run_Metadata"
                        },
                        {
                            "source": "#QualityFilterOuter.cwl/Split_Read_Pairs",
                            "id": "#QualityFilterOuter.cwl/Quality_Filter_Scatter/Split_Read_Pairs"
                        }
                    ]
                }
            ],
            "id": "#QualityFilterOuter.cwl"
        },
        {
            "inputs": [
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#SplitAndSubsample.cwl/Fastqs"
                },
                {
                    "type": {
                        "items": "string",
                        "type": "array"
                    },
                    "id": "#SplitAndSubsample.cwl/FilesToSkipSplitAndSubsample"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#SplitAndSubsample.cwl/NumRecordsPerSplit"
                },
                {
                    "type": "float",
                    "id": "#SplitAndSubsample.cwl/SubsampleRatio"
                },
                {
                    "type": "int",
                    "id": "#SplitAndSubsample.cwl/SubsampleSeed"
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "SplitAndSubsample splits, subsamples and formats read files to be deposited in QualityFilter.\n",
            "id": "#SplitAndSubsample.cwl",
            "steps": [
                {
                    "doc": "After scattering \"SplitAndSubsample\" on a File array, the output of each node is also an array. Thus, we are left with a nestled list. This JS expression flattens this list to deal with the split reads in PairReadFiles.cwl",
                    "out": [
                        "#SplitAndSubsample.cwl/FlattenOutput/SplitFastqList"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": {
                                    "items": {
                                        "items": "File",
                                        "type": "array"
                                    },
                                    "type": "array"
                                },
                                "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output/nestledSplitFastqList"
                            }
                        ],
                        "outputs": [
                            {
                                "type": {
                                    "items": "File",
                                    "type": "array"
                                },
                                "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output/SplitFastqList"
                            }
                        ],
                        "class": "ExpressionTool",
                        "expression": "${\n  return {SplitFastqList: [].concat.apply([], inputs.nestledSplitFastqList)}\n}\n",
                        "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output"
                    },
                    "id": "#SplitAndSubsample.cwl/FlattenOutput",
                    "in": [
                        {
                            "source": "#SplitAndSubsample.cwl/SplitAndSubsample/SplitAndSubsampledFastqs",
                            "id": "#SplitAndSubsample.cwl/FlattenOutput/nestledSplitFastqList"
                        }
                    ]
                },
                {
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "inputBinding": {
                                    "prefix": "--fastq-file-path"
                                },
                                "type": "File",
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/Fastq"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--files-to-skip-split-and-subsample",
                                    "itemSeparator": ","
                                },
                                "type": {
                                    "items": "string",
                                    "type": "array"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/FilesToSkipSplitAndSubsample"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--num-records"
                                },
                                "type": [
                                    "null",
                                    "long"
                                ],
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/NumRecordsPerSplit"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--subsample-ratio"
                                },
                                "type": "float",
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SubsampleRatio"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--subsample-seed"
                                },
                                "type": "int",
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SubsampleSeed"
                            }
                        ],
                        "requirements": [
                        ],
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "*.fastq.gz",
                                    "outputEval": "${ if (self.length === 0) { return [inputs.Fastq]; } else { return self; } }"
                                },
                                "type": {
                                    "items": "File",
                                    "type": "array"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SplitAndSubsampledFastqs"
                            },
                            {
                                "outputBinding": {
                                    "glob": "*.log"
                                },
                                "type": "File",
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/log"
                            }
                        ],
                        "baseCommand": [
                            "mist_split_fastq.py"
                        ],
                        "class": "CommandLineTool",
                        "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq"
                    },
                    "doc": "Allocate one docker/python process per file to do the actual file splitting.",
                    "scatter": [
                        "#SplitAndSubsample.cwl/SplitAndSubsample/Fastq"
                    ],
                    "in": [
                        {
                            "source": "#SplitAndSubsample.cwl/Fastqs",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/Fastq"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/FilesToSkipSplitAndSubsample",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/FilesToSkipSplitAndSubsample"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/NumRecordsPerSplit",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/NumRecordsPerSplit"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/SubsampleRatio",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/SubsampleRatio"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/SubsampleSeed",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/SubsampleSeed"
                        }
                    ],
                    "id": "#SplitAndSubsample.cwl/SplitAndSubsample",
                    "out": [
                        "#SplitAndSubsample.cwl/SplitAndSubsample/SplitAndSubsampledFastqs",
                        "#SplitAndSubsample.cwl/SplitAndSubsample/log"
                    ]
                }
            ],
            "outputs": [
                {
                    "outputSource": "#SplitAndSubsample.cwl/FlattenOutput/SplitFastqList",
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#SplitAndSubsample.cwl/SplitAndSubsampledFastqs"
                },
                {
                    "outputSource": "#SplitAndSubsample.cwl/SplitAndSubsample/log",
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#SplitAndSubsample.cwl/log"
                }
            ],
            "class": "Workflow"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#SubsampleSettings.cwl/_Subsample_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#SubsampleSettings.cwl/_Subsample_Seed"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#SubsampleSettings.cwl/Subsample_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#SubsampleSettings.cwl/Subsample_Seed"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  var subsamplingOutputs = {\n    Subsample_Reads: inputs._Subsample_Reads,\n    Subsample_Seed: inputs._Subsample_Seed\n  }\n  return subsamplingOutputs;\n}",
            "id": "#SubsampleSettings.cwl"
        },
        {
            "inputs": [
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#UncompressDatatables.cwl/Compressed_Data_Table"
                },
                {
                    "type": "File",
                    "id": "#UncompressDatatables.cwl/Compressed_Expression_Matrix"
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ],
            "outputs": [
                {
                    "outputSource": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompressed_File",
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#UncompressDatatables.cwl/Uncompressed_Data_Tables"
                },
                {
                    "outputSource": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompressed_File",
                    "type": "File",
                    "id": "#UncompressDatatables.cwl/Uncompressed_Expression_Matrix"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "id": "#UncompressDatatables.cwl/Uncompress_Datatable",
                    "out": [
                        "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompressed_File"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "inputBinding": {
                                    "position": 1
                                },
                                "type": "File",
                                "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner/Compressed_File"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "stdout": "$(inputs.Compressed_File.nameroot)",
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "$(inputs.Compressed_File.nameroot)"
                                },
                                "type": "File",
                                "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner/Uncompressed_File"
                            }
                        ],
                        "baseCommand": [
                            "gunzip"
                        ],
                        "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner",
                        "arguments": [
                            {
                                "position": 0,
                                "valueFrom": "-c"
                            }
                        ],
                        "class": "CommandLineTool",
                        "hints": [
                        ]
                    },
                    "scatter": [
                        "#UncompressDatatables.cwl/Uncompress_Datatable/Compressed_File"
                    ],
                    "in": [
                        {
                            "source": "#UncompressDatatables.cwl/Compressed_Data_Table",
                            "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Compressed_File"
                        }
                    ]
                },
                {
                    "out": [
                        "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompressed_File"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "inputBinding": {
                                    "position": 1
                                },
                                "type": "File",
                                "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner/Compressed_File"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "stdout": "$(inputs.Compressed_File.nameroot)",
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "$(inputs.Compressed_File.nameroot)"
                                },
                                "type": "File",
                                "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner/Uncompressed_File"
                            }
                        ],
                        "baseCommand": [
                            "gunzip"
                        ],
                        "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner",
                        "arguments": [
                            {
                                "position": 0,
                                "valueFrom": "-c"
                            }
                        ],
                        "class": "CommandLineTool",
                        "hints": [
                        ]
                    },
                    "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix",
                    "in": [
                        {
                            "source": "#UncompressDatatables.cwl/Compressed_Expression_Matrix",
                            "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Compressed_File"
                        }
                    ]
                }
            ],
            "id": "#UncompressDatatables.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/RSEC_Reads_Fastq"
                },
                {
                    "inputBinding": {
                        "position": 2
                    },
                    "type": "string",
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/Read_Limit"
                },
                {
                    "inputBinding": {
                        "position": 3
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/VDJ_Version"
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
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_pruned.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl/PyirCall"
                }
            ],
            "baseCommand": [
                "AssembleAndAnnotate.sh"
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs.cwl",
            "class": "CommandLineTool",
            "hints": [
                {
                    "coresMin": 1,
                    "ramMin": 3200,
                    "class": "ResourceRequirement"
                }
            ]
        },
        {
            "inputs": [
                {
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
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
            "outputs": [
                {
                    "outputSource": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/PyirCall",
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/igCalls"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "run": "#VDJ_Assemble_and_Annotate_Contigs.cwl",
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG",
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
                    "hints": [
                        {
                            "coresMin": "$(inputs.num_cores)",
                            "class": "ResourceRequirement"
                        }
                    ],
                    "scatter": [
                        "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/RSEC_Reads_Fastq"
                    ],
                    "out": [
                        "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl/VDJ_Assemble_and_Annotate_Contigs_IG/PyirCall"
                    ]
                }
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs_IG.cwl"
        },
        {
            "inputs": [
                {
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
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
            "outputs": [
                {
                    "outputSource": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/PyirCall",
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/tcrCalls"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "run": "#VDJ_Assemble_and_Annotate_Contigs.cwl",
                    "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR",
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
                    "hints": [
                        {
                            "coresMin": "$(inputs.num_cores)",
                            "class": "ResourceRequirement"
                        }
                    ],
                    "scatter": [
                        "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/RSEC_Reads_Fastq"
                    ],
                    "out": [
                        "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl/VDJ_Assemble_and_Annotate_Contigs_TCR/PyirCall"
                    ]
                }
            ],
            "id": "#VDJ_Assemble_and_Annotate_Contigs_TCR.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--seq-metrics"
                    },
                    "type": "File",
                    "id": "#VDJ_Compile_Results.cwl/Seq_Metrics"
                },
                {
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--cell-type-mapping-fp"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/cellTypeMapping"
                },
                {
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--ignore",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_Compile_Results.cwl/chainsToIgnore"
                },
                {
                    "inputBinding": {
                        "position": 8,
                        "prefix": "--e-value-for-j"
                    },
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/evalueJgene"
                },
                {
                    "inputBinding": {
                        "position": 7,
                        "prefix": "--e-value-for-v"
                    },
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/evalueVgene"
                },
                {
                    "inputBinding": {
                        "position": 5
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/igCalls"
                },
                {
                    "inputBinding": {
                        "position": 9,
                        "prefix": "--metadata-fp"
                    },
                    "type": "File",
                    "id": "#VDJ_Compile_Results.cwl/metadata"
                },
                {
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--putative-cells-json-fp"
                    },
                    "type": "File",
                    "id": "#VDJ_Compile_Results.cwl/putativeCells"
                },
                {
                    "inputBinding": {
                        "position": 6
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/tcrCalls"
                },
                {
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--vdj-version"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/vdjVersion"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "doc": "VDJ data per cell, with distribution based error correction",
                    "outputBinding": {
                        "glob": "*_VDJ_perCell.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/vdjCellsDatatable"
                },
                {
                    "doc": "VDJ data per cell, including non-putative cells, no error correction applied",
                    "outputBinding": {
                        "glob": "*_VDJ_perCell_uncorrected.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/vdjCellsDatatableUncorrected"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_Dominant_Contigs.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/vdjDominantContigs"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_metrics.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/vdjMetricsCsv"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_metrics.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/vdjMetricsJson"
                },
                {
                    "outputBinding": {
                        "glob": "*_DBEC_cutoff.png"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#VDJ_Compile_Results.cwl/vdjReadsPerCellByChainTypeFigure"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_Unfiltered_Contigs.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Compile_Results.cwl/vdjUnfilteredContigs"
                }
            ],
            "baseCommand": [
                "mist_vdj_compile_results.py"
            ],
            "id": "#VDJ_Compile_Results.cwl",
            "class": "CommandLineTool",
            "hints": [
                {
                    "ramMin": 32000,
                    "class": "ResourceRequirement"
                }
            ]
        },
        {
            "inputs": [
                {
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_GatherCalls.cwl/theCalls"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "VDJ_GatherCalls collect the outputs from the multi-processed VDJ step into one file.\n",
            "id": "#VDJ_GatherCalls.cwl",
            "steps": [
                {
                    "out": [
                        "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gatheredCalls"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": [
                                    {
                                        "items": [
                                            "null",
                                            "File"
                                        ],
                                        "type": "array"
                                    }
                                ],
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR/theCalls"
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
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "*_constant_region_called_pruned.csv.gz",
                                    "outputEval": "${\n  if (self.size == 0) {\n    throw(\"No outputs from PyIR detected in VDJ_GatherCalls!\");\n  } else {\n    return(self);\n  }\n}"
                                },
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR/gatheredCalls"
                            }
                        ],
                        "class": "CommandLineTool",
                        "arguments": [
                            {
                                "shellQuote": false,
                                "valueFrom": "${\n  if (!inputs.theCalls[0] ) {\n    return (\"echo \\\"No outputs from PyIR detected in VDJ_GatherCalls\\\"\")\n  }\n  var inputFiles = \"\"\n  if (!inputs.theCalls[0].path.split(\"_PrunePyIR\")[1]){\n    inputFiles = \"zcat\"\n    for (var i = 0; i < inputs.theCalls.length; i++) {\n      inputFiles += \" \" + inputs.theCalls[i].path\n    }\n    inputFiles += \" | \"\n  } else {\n    inputFiles = \"zcat \" + inputs.theCalls[0].path.split(\"VDJ\")[0] + \"*\" + inputs.theCalls[0].path.split(\"_PrunePyIR\")[1].split(\"_Number_\")[0] + \"_Number_*.csv.gz | \"\n  }\n  var outputFileName = \"\\\"gzip > \" + inputs.theCalls[0].nameroot.split(\"_Number_\")[0] + \"_constant_region_called_pruned.csv.gz\" + \"\\\"\"\n  var awkCommand =  \"awk \\'NR==1{F=$1;print | \" + outputFileName + \" } $1!=F { print | \" + outputFileName + \" }\\' \"\n  var outputCommand = inputFiles + awkCommand\n  return (outputCommand)\n}"
                            }
                        ],
                        "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR"
                    },
                    "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls",
                    "in": [
                        {
                            "source": "#VDJ_GatherCalls.cwl/theCalls",
                            "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/theCalls"
                        }
                    ]
                }
            ],
            "outputs": [
                {
                    "outputSource": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gatheredCalls",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_GatherCalls.cwl/gatheredCalls"
                }
            ],
            "class": "Workflow"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/Valid_Reads_Fastq"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/num_valid_reads"
                },
                {
                    "type": "string",
                    "id": "#VDJ_Preprocess_Reads.cwl/vdj_type"
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
                    "envDef": [
                        {
                            "envName": "CORES_ALLOCATED_PER_CWL_PROCESS",
                            "envValue": "8"
                        }
                    ],
                    "class": "EnvVarRequirement"
                }
            ],
            "outputs": [
                {
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/RSEC_Reads_Fastq",
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/RSEC_Reads_Fastq"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_cores",
                    "id": "#VDJ_Preprocess_Reads.cwl/num_cores"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "outputSource": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_splits",
                    "id": "#VDJ_Preprocess_Reads.cwl/num_splits"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "run": "#VDJ_RSEC_Reads.cwl",
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads/RSEC_Reads_Fastq"
                    ],
                    "requirements": [
                        {
                            "coresMin": 8,
                            "ramMin": "${ var est_ram = 0.0006 * parseInt(inputs.num_valid_reads) + 2000; var buffer = 1.25; est_ram *= buffer; if (est_ram < 2000) return 2000; if (est_ram > 370000) return 370000; return parseInt(est_ram); }",
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_RSEC_Reads",
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
                    ]
                },
                {
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Valid_Reads",
                        "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Trim_Report"
                    ],
                    "in": [
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/Valid_Reads_Fastq",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads/Valid_Reads_Fastq"
                        }
                    ],
                    "run": "#VDJ_Trim_Reads.cwl",
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_Trim_Reads",
                    "hints": [
                        {
                            "coresMin": 8,
                            "class": "ResourceRequirement"
                        }
                    ]
                },
                {
                    "out": [
                        "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_splits",
                        "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_cores"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits/num_valid_reads"
                            },
                            {
                                "type": "string",
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits/vdj_type"
                            }
                        ],
                        "outputs": [
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits/num_cores"
                            },
                            {
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits/num_splits"
                            }
                        ],
                        "class": "ExpressionTool",
                        "expression": "${\n  var ram_per_instance = 192 * 1024;\n  var num_cores = 96;\n  if (inputs.vdj_type == \"BCR\") {\n    ram_per_instance = 144 * 1024;\n    num_cores = 72;\n  }\n  var ram_per_split = 3200;\n  var num_splits_per_instance = parseInt(ram_per_instance / ram_per_split);\n  var num_splits = num_splits_per_instance;\n\n  var num_reads = parseInt(inputs.num_valid_reads);\n  if (num_reads != null) {\n    if (num_reads > 100000000)\n      num_splits = num_splits_per_instance * 2;\n      num_cores = num_cores * 2;\n  }\n\n  return ({\"num_splits\": num_splits, \"num_cores\": num_cores});\n}",
                        "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/determine_num_splits"
                    },
                    "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits",
                    "in": [
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/num_valid_reads",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/num_valid_reads"
                        },
                        {
                            "source": "#VDJ_Preprocess_Reads.cwl/vdj_type",
                            "id": "#VDJ_Preprocess_Reads.cwl/VDJ_num_splits/vdj_type"
                        }
                    ]
                }
            ],
            "id": "#VDJ_Preprocess_Reads.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--vdj-valid-reads",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_RSEC_Reads.cwl/Valid_Reads"
                },
                {
                    "inputBinding": {
                        "prefix": "--num-splits"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_RSEC_Reads.cwl/num_splits"
                }
            ],
            "requirements": [
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*RSEC_Reads_Fastq_*.tar.gz"
                    },
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_RSEC_Reads.cwl/RSEC_Reads_Fastq"
                }
            ],
            "baseCommand": "mist_vdj_rsec_reads.py",
            "class": "CommandLineTool",
            "id": "#VDJ_RSEC_Reads.cwl"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#VDJ_Settings.cwl/_VDJ_Version"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
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
            "class": "ExpressionTool",
            "expression": "${\n  var vdjVersion = null;\n  if (!inputs._VDJ_Version) {\n    vdjVersion = null;}\n  else {\n    var _VDJ_Version = inputs._VDJ_Version.toLowerCase();\n    if (_VDJ_Version === \"human\" || _VDJ_Version === \"hs\" || _VDJ_Version === \"human vdj - bcr and tcr\") {\n      vdjVersion = \"human\";\n    } else if (_VDJ_Version === \"humanbcr\" || _VDJ_Version === \"human vdj - bcr only\") {\n      vdjVersion = \"humanBCR\";\n    } else if (_VDJ_Version === \"humantcr\" || _VDJ_Version === \"human vdj - tcr only\") {\n      vdjVersion = \"humanTCR\";\n    } else if (_VDJ_Version === \"mouse\" || _VDJ_Version === \"mm\" || _VDJ_Version === \"mouse vdj - bcr and tcr\") {\n      vdjVersion = \"mouse\";\n    } else if (_VDJ_Version === \"mousebcr\" || _VDJ_Version === \"mouse vdj - bcr only\") {\n      vdjVersion = \"mouseBCR\";\n    } else if (_VDJ_Version === \"mousetcr\" || _VDJ_Version === \"mouse vdj - tcr only\") {\n      vdjVersion = \"mouseTCR\";\n    } else {\n      vdjVersion = inputs._VDJ_Version;\n    }\n  }\n\n  return ({\n  VDJ_Version: vdjVersion,\n  })\n}",
            "id": "#VDJ_Settings.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Trim_Reads.cwl/Valid_Reads_Fastq"
                }
            ],
            "requirements": [
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "cutadapt.log"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Trim_Reads.cwl/Trim_Report"
                },
                {
                    "outputBinding": {
                        "glob": "*vdjtxt.gz"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_Trim_Reads.cwl/Valid_Reads"
                }
            ],
            "baseCommand": "VDJ_Trim_Reads.sh",
            "class": "CommandLineTool",
            "id": "#VDJ_Trim_Reads.cwl"
        },
        {
            "inputs": [],
            "requirements": [
            ],
            "stdout": "output.txt",
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "output.txt",
                        "loadContents": true,
                        "outputEval": "$(self[0].contents)"
                    },
                    "type": "string",
                    "id": "#Version.cwl/version"
                }
            ],
            "baseCommand": [
                "mist_version.py"
            ],
            "id": "#Version.cwl",
            "class": "CommandLineTool"
        }
    ],
    "$namespaces": {
        "sbg": "https://sevenbridges.com#",
        "arv": "http://arvados.org/cwl#"
    }
}