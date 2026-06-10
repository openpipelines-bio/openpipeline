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
    // Check for correctness of mitochondrial and ribosomal gene detection arguments
    | map { id, state ->
      def new_state = [:]
      if (state.obs_name_mitochondrial_fraction && !state.var_name_mitochondrial_genes) {
        throw new RuntimeException("Using --obs_name_mitochondrial_fraction requires --var_name_mitochondrial_genes.")
      }
      if (!state.obs_name_mitochondrial_fraction && state.var_name_mitochondrial_genes) {
        new_state.obs_name_mitochondrial_fraction = "fraction_${state.var_name_mitochondrial_genes}"
      }
      if ((state.min_fraction_mito != null || state.max_fraction_mito != null) && !state.var_name_mitochondrial_genes) {
        throw new RuntimeException("Enabling --min_fraction_mito or --max_fraction_mito requires --var_name_mitochondrial_genes.")
      }
      if (state.var_gene_names && !state.var_name_mitochondrial_genes && !state.var_name_ribosomal_genes) {
        System.err.println("Warning: --var_gene_names is only required for mitochondrial gene detection and does nothing while \
                           not also setting --var_name_mitochondrial_genes or --var_name_ribosomal_genes.")
      }
      if (state.mitochondrial_gene_regex && !state.var_name_mitochondrial_genes && !state.var_name_ribosomal_genes) {
        System.err.println("Warning: --mitochondrial_gene_regex is only required for mitochondrial gene detection and does \
                           nothing while not also setting --var_name_mitochondrial_genes or --var_name_ribosomal_genes.")
      }
      if (state.obs_name_ribosomal_fraction && !state.var_name_ribosomal_genes) {
        throw new RuntimeException("Using --obs_name_ribosomal_fraction requires --var_name_ribosomal_genes.")
      }
      if (!state.obs_name_ribosomal_fraction && state.var_name_ribosomal_genes) {
        new_state.obs_name_ribosomal_fraction = "fraction_${state.var_name_ribosomal_genes}"
      }
      if ((state.min_fraction_ribo != null || state.max_fraction_ribo != null) && !state.var_name_ribosomal_genes) {
        throw new RuntimeException("Enabling --min_fraction_ribo or --max_fraction_ribo requires --var_name_ribosomal_genes.")
      }
      [id, state + new_state]
    }
    | qc.run(
      key: "qc_rna",
      fromState: { id, state ->
        // The rna singlesample processing allows for the optional filtering based on:
        // 1) fraction of mitochondrial genes per cell
        // 2) fraction of ribosomal genes per cell
        // 3) percentile of counts per cell
        // The behavior of the QC component must be tuned to the presence of these arguments.
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
          "modality": "rna",
          "layer": state.layer,
        ]
        
        if (state.min_percentile_counts || state.max_percentile_counts) {
        // If percentile-based filtering is enabled, total counts per cell must be calculated,
        // together with their log1p transform. The quantile filter operates on the
        // "log1p_total_counts" column, so both must be requested explicitly here rather than
        // relying on the qc component defaults.
          args += [
            "output_obs_total_counts_vars": "total_counts",
            "log1p_transform": true
          ]
        }
        if (state.var_name_mitochondrial_genes) {
          // Check if user has defined var columns to calculate metrics
          def new_var_qc_metrics = state.var_qc_metrics != null ? state.var_qc_metrics : []
          assert new_var_qc_metrics instanceof List
          // Add the mitochondrial genes var column to the columns to calculate statistics for if set.
          new_var_qc_metrics = ((new_var_qc_metrics as Set) + [state.var_name_mitochondrial_genes]) as List
          
          args += [
            "var_qc_metrics": new_var_qc_metrics,
            "obs_name_mitochondrial_fraction": state.obs_name_mitochondrial_fraction,
            "var_gene_names": state.var_gene_names,
            "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
            "mitochondrial_gene_regex": state.mitochondrial_gene_regex
          ]
        }

        if (state.var_name_ribosomal_genes) {
          // Check if user has defined var columns to calculate metrics
          def new_var_qc_metrics = state.var_qc_metrics != null ? state.var_qc_metrics : []
          assert new_var_qc_metrics instanceof List
          // Add the ribosomal genes var column to the columns to calculate statistics for if set.
          new_var_qc_metrics = ((new_var_qc_metrics as Set) + [state.var_name_ribosomal_genes]) as List
          
          args += [
            "var_qc_metrics": new_var_qc_metrics,
            "obs_name_ribosomal_fraction": state.obs_name_ribosomal_fraction,
            "var_gene_names": state.var_gene_names,
            "var_name_ribosomal_genes": state.var_name_ribosomal_genes,
            "ribosomal_gene_regex": state.ribosomal_gene_regex
          ]
        }

        return args
      },
      toState: ["input": "output"]
    )
    | delimit_fraction.run(
      key: "rna_filter_mitochondrial",
      runIf: {id, state -> state.var_name_mitochondrial_genes},
      fromState: {id, state -> 
      [
        "input": state.input,
        "obs_name_filter": "filter_mitochondrial",
        "min_fraction": state.min_fraction_mito,
        "max_fraction": state.max_fraction_mito,
        "obs_fraction_column": state.obs_name_mitochondrial_fraction,
      ]
      },
      toState: ["input": "output"]
    )
    | delimit_fraction.run(
      key: "rna_filter_ribosomal",
      runIf: {id, state -> state.var_name_ribosomal_genes},
      fromState: {id, state -> 
      [
        "input": state.input,
        "obs_name_filter": "filter_ribosomal",
        "min_fraction": state.min_fraction_ribo,
        "max_fraction": state.max_fraction_ribo,
        "obs_fraction_column": state.obs_name_ribosomal_fraction,
      ]
      },
      toState: ["input": "output"]
    )
    // cell filtering based on quantiles
    | filter_with_quantile.run(
      key: "rna_filter_by_percentile",
      runIf: { id, state -> state.min_percentile_counts || state.max_percentile_counts },
      fromState: { id, state ->
        [
          "input": state.input,
          "obs_min_quantile": state.min_percentile_counts,
          "obs_max_quantile": state.max_percentile_counts
        ]
      },
      args: [
          // Quantile filtering is always performed on the log-transformed total counts,
          // which are requested explicitly from the qc component above
          // (output_obs_total_counts_vars + log1p_transform).
          "obs_column": "log1p_total_counts",
          "obs_log1p_transform": false,
          "obs_name_filter": "filter_with_percentile"
      ],
      toState: ["input": "output"]
    )
    // cell filtering based on counts
    | filter_with_counts.run(
      key: "rna_filter_with_counts",
      fromState: { id, state ->
        [
          "input": state.input,
          "layer": state.layer,
          "obs_name_filter": "filter_with_counts",
          "var_name_filter": "filter_with_counts",
          "min_counts": state.min_counts,
          "max_counts": state.max_counts,
          "min_genes_per_cell": state.min_genes_per_cell,
          "max_genes_per_cell": state.max_genes_per_cell,
          "min_cells_per_gene": state.min_cells_per_gene,
        ]
      },
      toState: ["input": "output"]
    )
    // Doublet calling. Scrublet is run before do_filter so that it scores cells
    // on the full population. Removing low-quality cells first biases scrublet's
    // expected doublet rate (estimated from the number of cells given) and the
    // simulated-doublet neighborhood, lowering the doublet scores.
    | filter_with_scrublet.run(
      runIf: { id, state ->
        !state.skip_scrublet_doublet_detection
      },
      fromState: [
        "input": "input",
        "layer": "layer",
        "do_subset": "scrublet_do_subset",
        "scrublet_score_threshold": "scrublet_score_threshold",
        "expected_doublet_rate": "scrublet_expected_doublet_rate",
        "stdev_doublet_rate": "scrublet_stdev_doublet_rate",
        "n_neighbors": "scrublet_n_neighbors",
        "sim_doublet_ratio": "scrublet_sim_doublet_ratio",
        "min_counts": "scrublet_min_counts",
        "min_cells": "scrublet_min_cells",
        "min_gene_variablity_percent": "scrublet_min_gene_variability_percent",
        "num_pca_components": "scrublet_num_pca_components",
        "distance_metric": "scrublet_distance_metric",
        "allow_automatic_threshold_detection_fail": "scrublet_allow_automatic_threshold_detection_fail"
      ],
      args: [output_compression: "gzip"],
      toState: ["input": "output"]
    )
    | do_filter.run(
      key: "rna_do_filter",
      fromState: {id, state ->
        // do_filter does not need a layer argument because it filters all layers
        // from a modality.
        def stateMapping = [
          input: state.input,
          var_filter: ["filter_with_counts"],
          output: state.workflow_output,
          output_compression: "gzip"
        ]
        def obs_filter = ["filter_with_counts"]
        if (state.var_name_mitochondrial_genes) {
          obs_filter += ["filter_mitochondrial"]
        }
        if (state.var_name_ribosomal_genes) {
          obs_filter += ["filter_ribosomal"]
        }
        if (state.min_percentile_counts || state.max_percentile_counts) {
          obs_filter += ["filter_with_percentile"]
        }
        stateMapping += ["obs_filter": obs_filter]
        return stateMapping
      },
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}