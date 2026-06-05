workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // Calculate the qc metrics required for filtering, once, on the raw counts:
    //   - total_counts (.obs)      -> used to filter cells on --min_counts/--max_counts
    //   - num_nonzero_vars (.obs)  -> used to filter cells on --min_proteins_per_cell/--max_proteins_per_cell
    //   - num_nonzero_obs (.var)   -> used to filter proteins on --min_cells_per_protein
    | calculate_qc_metrics.run(
      key: "prot_calculate_qc_metrics",
      runIf: { id, state ->
        state.min_counts != null || state.max_counts != null ||
        state.min_proteins_per_cell != null || state.max_proteins_per_cell != null ||
        state.min_cells_per_protein != null
      },
      fromState: { id, state ->
        def filter_on_counts = state.min_counts != null || state.max_counts != null
        def filter_on_proteins_per_cell = state.min_proteins_per_cell != null || state.max_proteins_per_cell != null
        def filter_on_cells_per_protein = state.min_cells_per_protein != null
        [
          "input": state.input,
          "modality": "prot",
          "layer": state.layer,
          // top_n_vars is left unset (null) on purpose: an empty list triggers a
          // Viash parsing bug (int('') on the empty multiple-value separator).
          "top_n_vars": null,
          "output_obs_total_counts_vars": filter_on_counts ? "total_counts" : null,
          "output_obs_num_nonzero_vars": filter_on_proteins_per_cell ? "num_nonzero_vars" : null,
          "output_var_num_nonzero_obs": filter_on_cells_per_protein ? "num_nonzero_obs" : null,
          "output_var_total_counts_obs": null,
          "output_var_obs_mean": null,
          "output_var_pct_dropout": null,
        ]
      },
      toState: ["input": "output"]
    )
    // cell filtering on total counts per cell (.obs total_counts)
    | delimit_counts.run(
      key: "prot_filter_total_counts",
      runIf: { id, state -> state.min_counts != null || state.max_counts != null },
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "obs_count_column": ["total_counts"],
          "obs_name_filter": ["filter_total_counts"],
          "min_count": state.min_counts,
          "max_count": state.max_counts,
        ]
      },
      toState: ["input": "output"]
    )
    // cell filtering on the number of proteins per cell (.obs num_nonzero_vars)
    | delimit_counts.run(
      key: "prot_filter_proteins_per_cell",
      runIf: { id, state -> state.min_proteins_per_cell != null || state.max_proteins_per_cell != null },
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "obs_count_column": ["num_nonzero_vars"],
          "obs_name_filter": ["filter_proteins_per_cell"],
          "min_count": state.min_proteins_per_cell,
          "max_count": state.max_proteins_per_cell,
        ]
      },
      toState: ["input": "output"]
    )
    // protein filtering on the number of cells per protein (.var num_nonzero_obs)
    | delimit_counts.run(
      key: "prot_filter_cells_per_protein",
      runIf: { id, state -> state.min_cells_per_protein != null },
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "var_count_column": ["num_nonzero_obs"],
          "var_name_filter": ["filter_cells_per_protein"],
          "min_count": state.min_cells_per_protein,
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      key: "prot_do_filter",
      fromState : { id, state ->
        // do_filter does not need a layer argument because it filters all layers
        // from a modality. Only the filter columns that were actually created
        // (i.e. for which a threshold was requested) are referenced here.
        def obs_filter = []
        if (state.min_counts != null || state.max_counts != null) {
          obs_filter += ["filter_total_counts"]
        }
        if (state.min_proteins_per_cell != null || state.max_proteins_per_cell != null) {
          obs_filter += ["filter_proteins_per_cell"]
        }
        def var_filter = []
        if (state.min_cells_per_protein != null) {
          var_filter += ["filter_cells_per_protein"]
        }
        def newState = [
          "input": state.input,
          "modality": "prot",
          "output_compression": "gzip",
          "output": state.workflow_output
        ]
        // Only pass the filter arguments when columns were actually created;
        // an empty list would be passed as a single empty string to do_filter.
        if (obs_filter) {
          newState["obs_filter"] = obs_filter
        }
        if (var_filter) {
          newState["var_filter"] = var_filter
        }
        return newState
      },
      toState: ["output": "output"],
    )
    | setState(["output"])

  emit:
  output_ch
}
