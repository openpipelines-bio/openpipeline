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
          "var_qc_metrics": null
        ]
        newState
      }
    )
    | setState(["output"])


  emit:
  output_ch
}