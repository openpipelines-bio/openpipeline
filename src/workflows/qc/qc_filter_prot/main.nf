workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Calculate the total counts (and its log1p) per cell, used by the count and
    // quantile filters below. Other QC metrics are disabled.
    | calculate_qc_metrics.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "layer": state.layer,
          "top_n_vars": null,
          "output_obs_num_nonzero_vars": null,
          "output_obs_total_counts_vars": "total_counts",
          "output_var_num_nonzero_obs": null,
          "output_var_total_counts_obs": null,
          "output_var_obs_mean": null,
          "output_var_pct_dropout": null,
          "log1p_transform": true,
        ]
      },
      toState: ["input": "output"]
    )
    // Flag cells with too few total counts.
    | filter_with_counts.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "layer": state.layer,
          "obs_name_filter": "filter_counts_prot",
          "min_counts": state.min_count,
        ]
      },
      toState: ["input": "output"]
    )
    // Flag cells above the upper quantile of (log1p) total counts.
    | filter_with_quantile.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "obs_column": "total_counts",
          "obs_log1p_transform": true,
          "obs_max_quantile": state.max_quantile,
          "obs_name_filter": "filter_quantile_prot",
        ]
      },
      toState: ["input": "output"]
    )
    // Emit the flagged (but not subsetted) data. Subsetting and cross-modality
    // intersection are handled by the top-level qc_filter workflow.
    | setState(["output": "input"])

  emit:
  output_ch
}
