workflow run_wf {
  take:
  input_ch

  main:

  output_ch = input_ch
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | clr.run(
      fromState: [
        "input": "input",
        "input_layer": "layer",
        ],
      toState: ["input": "output"],
      args: [ 
        output_layer: "clr", 
        modality: "prot"
      ]
    )
    | prot_qc.run(
      // TODO: remove when viash 0.8.3 is released
      key: "prot_qc",
      fromState: { id, state ->
        def newState = [
          "id": id,
          "output": state.workflow_output,
          "input": state.input,
          "top_n_vars": state.top_n_vars,
          "var_qc_metrics": null,
          "input_layer": state.layer, // Use the non-transformed layer
          "modality": "prot",
          "var_name_mitochondrial_genes": null,
          "output_obs_num_nonzero_vars": state.output_obs_num_nonzero_vars,
          "output_obs_total_counts_vars": state.output_obs_total_counts_vars,
          "num_nonzeoutput_var_num_nonzero_obsro_obs": state.output_var_num_nonzero_obs,
          "output_var_total_counts_obs": state.output_var_total_counts_obs,
          "output_var_obs_mean": state.output_var_obs_mean,
          "output_var_pct_dropout": state.pct_dropout
        ]
        newState
      }
    )
    | setState(["output"])


  emit:
  output_ch
}