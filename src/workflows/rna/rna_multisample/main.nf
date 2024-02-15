workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | normalize_total.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "input_layer": state.layer, 
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
    | highly_variable_features_scanpy.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "layer": "log_normalized",
          "modality": state.modality,
          "var_name_filter": state.highly_variable_features_var_output,
          "n_top_features": state.highly_variable_features_n_top_features,
          "flavor": state.highly_variable_features_flavor,
          "obs_batch_key": state.highly_variable_features_obs_batch_key
        ]
      },
      toState: ["input": "output"],
    )
    | rna_qc.run(
      // TODO: remove when viash 0.8.3 is released
      key: "rna_qc",
      fromState: {id, state ->
        [
          "id": id,
          "input": state.input,
          "output": state.workflow_output,
          "layer": state.layer, // Use the non-transformed layer
          "output_compression": "gzip",
          "modality": state.modality,
          "var_qc_metrics": state.var_qc_metrics,
          "top_n_vars": state.top_n_vars,
          "output_obs_num_nonzero_vars": state.output_obs_num_nonzero_vars,
          "output_obs_total_counts_vars": state.output_obs_total_counts_vars,
          "output_var_num_nonzero_obs": state.output_var_num_nonzero_obs,
          "output_var_total_counts_obs": state.output_var_total_counts_obs,
          "output_var_obs_mean": state.output_var_obs_mean,
          "output_var_pct_dropout": state.output_var_pct_dropout
        ]
      },
    )
    | setState(["output"])

  emit:
  output_ch
}