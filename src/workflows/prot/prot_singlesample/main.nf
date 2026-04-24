workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
  | qc.run(
    key: "qc_prot",
    fromState: { id, state ->
      // The prot singlesample processing allows for optional filtering based on the percentile distribution of counts
      // The behavior of the QC component must be tuned to the presence of the percentile-based filtering arguments as well
      def args = [
        "id": id,
        "input": state.input,
        // disable all qc metrics by default
        "top_n_vars": [],
        "output_obs_num_nonzero_vars": null,
        "output_obs_total_counts_vars": null,
        "output_var_num_nonzero_obs": null,
        "output_var_total_counts_obs": null,
        "output_var_obs_mean": null,
        "output_var_pct_dropout": null,
        "output": state.output,
        "modality": "prot",
        "layer": state.layer,
      ]
        
      if (state.min_percentile_counts || state.max_percentile_counts) {
      // If percentile-based filtering is enabled, total counts per cell must be calculated.
        args += [
          "output_obs_total_counts_vars": "total_counts"
        ]
      }

      return args
      },
      toState: ["input": "output"]
    )
    // filtering
    | filter_with_quantile.run(
      key: "prot_filter_by_percentile",
      runIf: { id, state -> state.min_percentile_counts || state.max_percentile_counts },
      fromState: { id, state ->
        [
          "input": state.input,
          "obs_min_quantile": state.min_percentile_counts,
          "obs_max_quantile": state.max_percentile_counts,
          "obs_log1p_transform": state.log_transform_total_counts
        ]
      },
      args: [
          "modality": "prot",
          "obs_column": "total_counts",
          "obs_name_filter": "filter_with_percentile"
      ],
      toState: ["input": "output"]
    )
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
        def obs_filter = ["filter_with_counts"]
        if (state.min_percentile_counts || state.max_percentile_counts) {
          obs_filter += ["filter_with_percentile"]
        }
        def newState = [
          "input": state.input,
          "obs_filter": obs_filter,
          "modality": "prot",
          "var_filter": "filter_with_counts",
          "output_compression": "gzip",
          "output": state.workflow_output
        ]
        return newState
      },
      toState: ["output": "output"],
    )
    | setState(["output"])

  emit:
  output_ch
}
