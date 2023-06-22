nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { normalize_total } from targetDir + '/transform/normalize_total/main.nf'
include { log1p } from targetDir + '/transform/log1p/main.nf'
include { filter_with_hvg } from targetDir + '/filter/filter_with_hvg/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'
include { calculate_qc_metrics } from targetDir + '/qc/calculate_qc_metrics/main.nf'
include { delete_layer } from targetDir + '/transform/delete_layer/main.nf'
include { add_id } from targetDir + "/metadata/add_id/main.nf"

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
  parsed_arguments_ch = input_ch
    | preprocessInputs("config": config)
    // add the id to the arguments
    | pmap { id, data ->
      def new_data = data + [ input_id: data.sample_id ]
      [id, new_data]
    }
    | setWorkflowArguments (
      "add_id": ["input": "input",
                 "input_id": "sample_id",
                 "obs_output": "add_id_obs_output",
                 "make_observation_keys_unique": "add_id_make_observation_keys_unique"],
      "concat": ["input_id": "sample_id"],
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
    | getWorkflowArguments(key: "add_id")
    
  add_id_ch = parsed_arguments_ch
    | filter{ it[1].add_id_to_obs }
    // The add_id processes will be executed several times
    // Before that, we must make the ID (first element of list)
    // unique. The global ID is stored so that it can be retreived later.
    | pmap {id, data, other ->
      [id, data, other + [id: id]]
    }
    // Split the input into multiple events in the channel so that 
    // the add_id process can be run multiple times
    | pFlatMap { id, data ->
      def singleInputs = [data.input_id, data.input]
        .transpose()
        .collect({list -> ["input_id": list[0], "input": list[1]]})
      def result = singleInputs.collect({map -> 
        [id + "_" + map.input_id, // Make the IDs unique
         map + data.findAll{!(['input_id', 'input'].contains(it.key))},
        ]
      })
      result
    }
    | add_id.run(auto: [ simplifyOutput: false ])
    | collect(sort: true, flat: false) // Join the results of the different events
    // Reformat the event to a proper single event again
    // and add the original ID
    | map { list -> 
      def other_arguments = list[0][2]
      def passthrough = list[0].drop(3)
      def globalID = other_arguments["id"]
      def inputs = list.collect({it -> it[1].output})
      [globalID, ["input": inputs], other_arguments] + passthrough
     }

    // Do nothing for the samples that do not need to have their ID 
    // added to the MuData object
  no_id_added_ch = parsed_arguments_ch
    | filter{ ! (it[1].add_id_to_obs) }

  samples_with_id_ch = add_id_ch.mix(no_id_added_ch)

  output_ch = samples_with_id_ch
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
      args: [ layer: "log_normalized", output_compression: "gzip"]
    )
    | getWorkflowArguments(key: "qc_metrics")
    | calculate_qc_metrics.run(
      // layer: null to use .X and not log transformed
      args: [
        input_layer: null,
        output_compression: "gzip"
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
    add_id_to_obs: false
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

workflow test2_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "combined_samples_rna_add_id",
    sample_id: "mouse;brain",
    input: params.resources_test + "/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu;" + params.resources_test + "/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
    add_id_make_observation_keys_unique: true,
    add_id_to_obs: true,
    add_id_obs_output: "foo_column"
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
      assert output_list[0][0] == "combined_samples_rna_add_id" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}