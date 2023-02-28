nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { normalize_total } from targetDir + '/transform/normalize_total/main.nf'
include { log1p } from targetDir + '/transform/log1p/main.nf'
include { filter_with_hvg } from targetDir + '/filter/filter_with_hvg/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'
include { calculate_qc_metrics } from targetDir + '/qc/calculate_qc_metrics/main.nf'
include { delete_layer } from targetDir + '/transform/delete_layer/main.nf'

include { readConfig; helpMessage; channelFromParams; preprocessInputs } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/rna_multisample/config.vsh.yaml")

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
    // add the id to the arguments
    | pmap { id, data ->
      def new_data = data + [ input_id: data.sample_id ]
      [id, new_data]
    }
    | setWorkflowArguments (
      "concat": [:],
      "normalize_total": [:],
      "log1p": [:],
      "delete_layer": [:],
      "filter_with_hvg": [
        "var_name_filter": "filter_with_hvg_var_output",
        "n_top_genes": "filter_with_hvg_n_top_genes",
        "flavor": "filter_with_hvg_flavor",
        "obs_batch_key": "filter_with_hvg_obs_batch_key"
      ],
      "qc_metrics": [
        "var_qc_metrics": "var_qc_metrics",
        "top_n_vars": "top_n_vars"
      ]
    )
    | getWorkflowArguments(key: "concat")
    | concat
    | getWorkflowArguments(key: "normalize_total")
    | normalize_total.run( 
      args: [ output_layer: "normalized" ]
    )

    | getWorkflowArguments(key: "log1p")
    | log1p.run( 
      args: [ output_layer: "log_normalized", input_layer: "normalized" ]
    )

    | getWorkflowArguments(key: "delete_layer")
    | delete_layer.run(
      args: [ layer: "normalized", modality: "rna" ]
    )

    // feature annotation
    | getWorkflowArguments(key: "filter_with_hvg")
    | filter_with_hvg.run(
      auto: [ publish: true ],
      args: [ layer: "log_normalized"]
    )
    | getWorkflowArguments(key: "qc_metrics")
    | calculate_qc_metrics.run(
      // layer: null to use .X and not log transformed
      args: [
        input_layer: null,       
      ]
    )
    | map {list -> [list[0], list[1]] + list.drop(3)}

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "combined_samples_rna",
    sample_id: "mouse;brain",
    input: params.resources_test + "/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu;" + params.resources_test + "/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
  ]

  output_ch =
    channelFromParams(testParams, config)
    | map {list -> list + [test_passthrough: "test"]}
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain three elements; [id, file, passthrough]"
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