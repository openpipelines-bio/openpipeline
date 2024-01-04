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
      fromState: ["input": "input"],
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
          "input_layer": null, // layer: null to use .X and not log transformed
          "modality": "prot",
          "var_name_mitochondrial_genes": null,
          "num_nonzero_vars": state.num_nonzero_vars,
          "total_counts_var": state.total_counts_var,
          "num_nonzero_obs": state.num_nonzero_obs,
          "total_counts_obs": state.total_counts_obs,
          "obs_mean": state.obs_mean,
          "pct_dropout": state.pct_dropout
        ]
        newState
      }
    )
    | setState(["output"])


  emit:
  output_ch
}