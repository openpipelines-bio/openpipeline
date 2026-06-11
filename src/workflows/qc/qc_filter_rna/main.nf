workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Detect mitochondrial genes and compute the per-cell mitochondrial fraction.
    | grep_annotation_column.run(
      key: "grep_mitochondrial_genes",
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "input_layer": state.layer,
          "input_column": state.var_gene_names,
          "matrix": "var",
          "regex_pattern": state.mitochondrial_gene_regex,
          "output_match_column": "mitochondrial",
          "output_fraction_column": "fraction_mitochondrial",
        ]
      },
      toState: ["input": "output"]
    )
    // Calculate the total counts (and its log1p) per cell, used by the count and
    // quantile filters below. Other QC metrics are disabled.
    | calculate_qc_metrics.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
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
    // Flag cells with too few total counts (empty droplets / low quality).
    | filter_with_counts.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "layer": state.layer,
          "obs_name_filter": "filter_counts_rna",
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
          "modality": "rna",
          "obs_column": "total_counts",
          "obs_log1p_transform": true,
          "obs_max_quantile": state.max_quantile,
          "obs_name_filter": "filter_quantile_rna",
        ]
      },
      toState: ["input": "output"]
    )
    // Flag cells with a mitochondrial fraction above the threshold.
    | delimit_fraction.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "obs_fraction_column": "fraction_mitochondrial",
          "max_fraction": state.max_fraction_mito,
          "obs_name_filter": "filter_mito_rna",
        ]
      },
      toState: ["input": "output"]
    )
    // Doublet detection (flag only). Skipped when --skip_scrublet is set.
    | filter_with_scrublet.run(
      runIf: {id, state -> !state.skip_scrublet},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "layer": state.layer,
          "obs_name_filter": "filter_scrublet",
          "scrublet_score_threshold": state.scrublet_score_threshold,
          "expected_doublet_rate": state.scrublet_expected_doublet_rate,
          "min_counts": state.scrublet_min_counts,
          "min_cells": state.scrublet_min_cells,
          "min_gene_variablity_percent": state.scrublet_min_gene_variability_percent,
          "num_pca_components": state.scrublet_num_pca_components,
        ]
      },
      args: [output_compression: "gzip"],
      toState: ["input": "output"]
    )
    // Emit the flagged (but not subsetted) data. Subsetting and cross-modality
    // intersection are handled by the top-level qc_filter workflow.
    | setState(["output": "input"])

  emit:
  output_ch
}
