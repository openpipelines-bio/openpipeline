workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // filtering
    | filter_with_counts.run(
      key: "prot_filter_with_counts",
      fromState: { id, state ->
        def newState = [
          "input": state.input,
          "min_counts": state.min_counts,
          "max_counts": state.max_counts,
          "min_genes_per_cell": state.min_proteins_per_cell,
          "max_genes_per_cell": state.max_proteins_per_cell,
          "min_cells_per_gene": state.min_cells_per_protein,
          "obs_name_filter": "filter_with_counts",
          "var_name_filter": "filter_with_counts",
          "modality": "prot",
          "layer": state.layer,
        ]
        newState
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      key: "prot_do_filter",
      fromState : { id, state ->
        // do_filter does not need a layer argument because it filters all layers
        // from a modality.
        def newState = [
          "input": state.input,
          "obs_filter": "filter_with_counts",
          "modality": "prot",
          "var_filter": "filter_with_counts",
          "output_compression": "gzip",
          "output": state.workflow_output
        ]
        return newState
      },
      toState: ["output": "output"],
      auto: [ publish: true ]
    )
    | setState(["output"])

  emit:
  output_ch
}
