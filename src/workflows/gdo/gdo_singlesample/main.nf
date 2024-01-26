workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Set aside the output for this workflow to avoid conflicts
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // cell filtering
    | filter_with_counts.run(
      key: "gdo_filter_with_counts",
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": "gdo",
          "obs_name_filter": "filter_with_counts",
          "var_name_filter": "filter_with_counts",
          "min_counts": state.min_counts,
          "max_counts": state.max_counts,
          "min_genes_per_cell": state.min_guides_per_cell,
          "max_genes_per_cell": state.max_guides_per_cell,
          "min_cells_per_gene": state.min_cells_per_guide,
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      key: "gdo_do_filter",
      fromState : { id, state ->
        def newState = [
          "input": state.input,
          "obs_filter": "filter_with_counts",
          "modality": "gdo",
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