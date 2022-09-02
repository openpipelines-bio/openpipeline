nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { normalize_total } from targetDir + '/transform/normalize_total/main.nf'
include { log1p } from targetDir + '/transform/log1p/main.nf'
include { filter_with_hvg } from targetDir + '/filter/filter_with_hvg/main.nf'
include { concat } from targetDir + '/integrate/concat/main.nf'
include { delete_layer } from targetDir + '/transform/delete_layer/main.nf'

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/process_rna/multisample/config.vsh.yaml")

// keep track of whether this is an integration test or not
global_params = [ do_publish: true ]

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, list -> ["combined_samples_rna", 
                          ["input": list.get("input"),
                           "sample_names":  list.get("id")]
                    ]}
    | concat
    // normalisation
    | normalize_total.run( 
      args: [ output_layer: "normalized" ]
    )
    | log1p.run( 
      args: [ output_layer: "log_normalized", input_layer: "normalized" ]
    )
    | delete_layer.run(
      args: [ layer: "normalized", modality: "rna" ]
    )
    // feature annotation
    | filter_with_hvg.run(
      auto: [ publish: global_params.do_publish ],
      args: [ layer: "log_normalized"]
    )

  emit:
  output_ch
}

workflow test_wf {
  // don't publish output
  global_params.do_publish = false

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "mouse;brain",
    input: params.resources_test + "/concat/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu;" + params.resources_test + "/concat/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu",
  ]

  output_ch =
    viashChannel(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain three elements; [id, file]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "combined_samples_rna" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}