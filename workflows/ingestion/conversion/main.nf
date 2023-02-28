nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf" 
include { from_10xmtx_to_h5mu } from targetDir + "/convert/from_10xmtx_to_h5mu/main.nf" 
include { from_h5ad_to_h5mu } from targetDir + "/convert/from_h5ad_to_h5mu/main.nf" 
include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"


config = readConfig("$workflowDir/ingestion/conversion/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
}

workflow run_wf {
    take:
        input_ch

    main:
        preprocessed_ch = input_ch
          | preprocessInputs("config": config)
        commonOptions = [
            auto: [ publish: true ]
        ]
        ch_10xh5 = preprocessed_ch
            | filter{ it[1].input_type == "10xh5" }
            | from_10xh5_to_h5mu.run(commonOptions)

        ch_10xmtx = preprocessed_ch
            | filter{ it[1].input_type  == "10xmtx" }
            | from_10xmtx_to_h5mu.run(commonOptions)

        ch_h5ad = preprocessed_ch
            | filter{ it[1].input_type  == "h5ad" }
            | from_h5ad_to_h5mu.run(commonOptions)

        
        // /* Combine the different conversion channels */
        all_ch = ch_10xh5 
            | mix( ch_10xmtx, ch_h5ad )
        output_ch = all_ch
    emit:
        output_ch
}

workflow test_wf {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"

  testParams = [
    param_list: [
      [
        id: "10xh5_test",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
        input_type: "10xh5",
        modality: null
      ],
      [
        id: "10xmtx_test",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix",
        input_type: "10xmtx",
        modality: null,
        output: "\$id.h5mu"
      ],
      [
        id: "h5ad",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix",
        input_type: "10xmtx",
        modality: "rna",
        output: "\$key.h5mu"
      ]
    ]
  ]
    
  output_ch =
    channelFromParams(testParams, config)
    | view { "Input: $it" }    
    | run_wf
    | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        assert output_list.size() == 3 : "output channel should contain three events"
      }
}