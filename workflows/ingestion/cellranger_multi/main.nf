nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"


include { cellranger_multi } from targetDir + "/mapping/cellranger_multi/main.nf"
include { from_cellranger_multi_to_h5mu } from targetDir + "/convert/from_cellranger_multi_to_h5mu/main.nf"

include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/ingestion/cellranger_multi/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs("config": config)
    // split params for downstream components
    | setWorkflowArguments(
      cellranger_multi: [
        "input": "input",
        "cell_multiplex_sample_id": "cell_multiplex_sample_id",
        "cell_multiplex_oligo_ids": "cell_multiplex_oligo_ids",
        "cell_multiplex_description": "cell_multiplex_description",
        "gex_expect_cells": "gex_expect_cells",
        "gex_chemistry": "gex_chemistry",
        "gex_secondary_analysis": "gex_secondary_analysis",
        "gex_generate_bam": "gex_generate_bam",
        "gex_include_introns": "gex_include_introns",
        "library_id": "library_id",
        "library_type": "library_type",
        "library_subsample": "library_subsample",
        "library_lanes": "library_lanes"
      ],
      from_cellranger_multi_to_h5mu: [
        "output": "output_h5mu",
        "uns_metrics": "uns_metrics"
      ]
    )

    | getWorkflowArguments(key: "cellranger_multi")
    | cellranger_multi.run(auto: [ publish: true ])
    | from_cellranger_multi_to_h5mu.run(auto: [ publish: true ])
     | pmap { id, data, output_data ->
      [ id, output_data + [h5mu: data] ]
    }


  emit:
  output_ch
}


workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input:[params.resources_test + "/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R1_001.fastq.gz",
           params.resources_test + "/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R2_001.fastq.gz",
           params.resources_test + "/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R1_001.fastq.gz",
           params.resources_test + "/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R2_001.fastq.gz",
           params.resources_test + "/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R1_001.fastq.gz",
           params.resources_test + "/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R2_001.fastq.gz"],
    gex_reference: params.resources_test + "/reference_gencodev41_chr1/reference_cellranger.tar.gz",
    vdj_reference: params.resources_test + "/10x_5k_anticmv/raw/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz",
    feature_reference: params.resources_test + "/10x_5k_anticmv/raw/feature_reference.csv",
    library_id: ["5k_human_antiCMV_T_TBNK_connect_GEX_1_subset", "5k_human_antiCMV_T_TBNK_connect_AB_subset", "5k_human_antiCMV_T_TBNK_connect_VDJ_subset"],
    library_type: ["Gene Expression", "Antibody Capture", "VDJ"]
  ]

  output_ch =
    channelFromParams(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}