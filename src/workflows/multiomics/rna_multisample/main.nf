workflow run_wf {
  take:
  input_ch

  main:
  parsed_arguments_ch = input_ch
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }

  // Check that the output name is unique for all samples (we concat into one file)
  parsed_arguments_ch
    | toSortedList
    | map { list ->
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