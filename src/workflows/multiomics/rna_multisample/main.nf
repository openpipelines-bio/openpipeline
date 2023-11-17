nextflow.enable.dsl=2

workflowDir = params.rootDir + "/src/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { normalize_total } from targetDir + '/transform/normalize_total/main.nf'
include { log1p } from targetDir + '/transform/log1p/main.nf'
include { filter_with_hvg } from targetDir + '/filter/filter_with_hvg/main.nf'
include { calculate_qc_metrics } from targetDir + '/qc/calculate_qc_metrics/main.nf'
include { delete_layer } from targetDir + '/transform/delete_layer/main.nf'
include { add_id } from targetDir + "/metadata/add_id/main.nf"

include { readConfig; helpMessage; channelFromParams; preprocessInputs } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/rna_multisample/config.vsh.yaml")

workflow rna_multisample {
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
    | pmap {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }

  // Check that the output name is unique for all samples (we concat into one file)
  parsed_arguments_ch
    | toSortedList
    | map  { list ->
        found_output_files = list.collect{it[1].getOrDefault('workflow_output', null)}.unique()
        assert found_output_files.size() < 2, "The specified output file is not the same for all samples. Found: $found_output_files"
    }

  output_ch = parsed_arguments_ch
    | normalize_total.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "output_layer": "normalized",
          "modality": state.modality
        ]
      },
      toState: ["input": "output"],
    )
    | log1p.run( 
      fromState: { id, state ->
        [
          "input": state.input,
          "output_layer": "log_normalized",
          "input_layer": "normalized",
          "modality": state.modality
        ]
      },
      toState: ["input": "output"]
    )
    | delete_layer.run(
      fromState: {id, state -> 
        [
          "input": state.input,
          "layer": "normalized",
          "modality": state.modality
        ]
      },
      toState: ["input": "output"]
    )
    | filter_with_hvg.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "layer": "log_normalized",
          "modality": state.modality,
          "var_name_filter": state.filter_with_hvg_var_output,
          "n_top_genes": state.filter_with_hvg_n_top_genes,
          "flavor": state.filter_with_hvg_flavor,
          "obs_batch_key": state.filter_with_hvg_obs_batch_key
        ]
      },
      toState: ["input": "output"],
    )
    | calculate_qc_metrics.run(
      // layer: null to use .X and not log transformed
      fromState: {id, state ->
        [
          "input": state.input,
          "output": state.workflow_output,
          "input_layer": null,
          "output_compression": "gzip",
          "modality": state.modality,
          "var_qc_metrics": state.var_qc_metrics,
          "top_n_vars": state.top_n_vars,
        ]
      },
      toState: {id, output, state -> 
        ["output": output.output]
      },
      key: "rna_calculate_qc_metrics",
      auto: [ publish: true ]
    )

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
      [id: "concatenated_file",
       input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
       output: "concatenated_file.final.h5mu"
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    | map {list -> list + [test_passthrough: "test"]}
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain three elements; [id, file, passthrough]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "concatenated_file" : "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["concatenated_file.final.h5mu"] as Set)
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}
