nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { publish }  from targetDir + '/transfer/publish/main.nf'
include { grep_annotation_column }  from targetDir + '/metadata/grep_annotation_column/main.nf'
include { calculate_qc_metrics } from targetDir + '/qc/calculate_qc_metrics/main.nf'


include { readConfig; helpMessage; readCsv; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { strictMap as smap; passthroughMap as pmap; } from workflowDir + "/utils/DataflowHelper.nf"
config = readConfig("$workflowDir/qc/qc/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf

}

workflow qc {
  take:
  input_ch

  main:
    preproc_ch = input_ch
    // Avoid conflict between output from component and output for this workflow
    | pmap {id, state -> 
      assert state.output, "Output must be defined"
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // Add default for var_qc_metrics component
    | pmap {id, state ->
      def var_qc_default = []
      // Remove the var_qc_metric argument from the state if its value is null (not specified)
      def new_state = state.findAll { it.key != "var_qc_metrics" || it.value == null }
      if (state.var_name_mitochondrial_genes) {
        var_qc_default.add(state.var_name_mitochondrial_genes)
      }
      // Get the new state, but make sure to overwrite var_qc_metrics if the user has set it.
      new_state = ["var_qc_metrics": var_qc_default.join(",")] + new_state
      [id, new_state]
    }

    with_grep_ch = preproc_ch
    | filter { it -> it[1].var_name_mitochondrial_genes }
    | grep_annotation_column.run(
      fromState: { id, state ->
        def stateMapping = [
          "input": state.input,
          "modality": state.modality,
          "input_column": state.var_gene_names,
          "matrix": "var",
          "output_match_column": state.var_name_mitochondrial_genes,
          "regex_pattern": state.mitochondrial_gene_regex
        ]
        stateMapping.output_fraction_column = state.obs_name_mitochondrial_fraction ? state.obs_name_mitochondrial_fraction: "fraction_$state.var_name_mitochondrial_genes"
        return stateMapping
      },
      toState: ["input": "output"]
    )

    without_grep_ch = preproc_ch
      | filter { it -> !it[1].var_name_mitochondrial_genes }

    output_ch = without_grep_ch.mix(with_grep_ch)
      | calculate_qc_metrics.run(
        fromState: { id, state ->
          def newState = [
            "input": state.input,
            "modality": state.modality,
            "layer": state.layer,
            "top_n_vars": state.top_n_vars,
            "var_qc_metrics_fill_na_value": state.var_qc_metrics_fill_na_value
          ]
          if (state.var_qc_metrics) {
            newState += ["var_qc_metrics": state.var_qc_metrics]
          }
        return newState
        },
        // use map when viash 0.7.6 is released
        // related to https://github.com/viash-io/viash/pull/515
        toState: ["input": "output"]
      )
      | publish.run(
        fromState: { id, state -> [
            "input": state.input,
            "output": state.workflow_output,
            "compression": "gzip"
          ]
        },
        auto: [ publish: true ]
      )      

  emit:
  output_ch
}

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    | preprocessInputs("config": config)
    | qc

  emit:
    output_ch
}


// ===============================
// === start of test workflows ===
// ===============================

workflow test_wf {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
        [
          id: "mouse",
          input: params.resources_test + "/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
          publish_dir: "foo/",
          output: "qc_metrics_mouse.h5mu"
        ],
        [
          id: "human",
          input: params.resources_test + "/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
          publish_dir: "foo/",
          output: "qc_metrics_human.h5mu"
        ]
      ],
      var_name_mitochondrial_genes: "mitochondrial"
    ]


  output_ch =
    channelFromParams(testParams, config)
      | view { "Input: $it" }
      | run_wf
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        assert output_list.size() == 2 : "output channel should contain two events"
        assert (output_list.collect({it[0]}) as Set).equals(["mouse", "human"] as Set): "Output ID should be same as input ID"
        assert (output_list.collect({it[1].output.getFileName().toString()}) as Set) == ["qc_metrics_mouse.h5mu", "qc_metrics_human.h5mu"] as Set: "Output files should be named qc_metrics_mouse.h5mu and qc_metrics_human.h5mu"
      }
  
}
