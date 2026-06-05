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
      fromState: { id, state ->
        // The rna singlesample processing allows detecting mitochondrial genes and filtering based
        // on the fraction of mitochondrial genes per cell
        // This behaviour is optional based on the presence of var_name_mitochondrial_genes
        // The behavior of other components must be tuned to this argument as well
        // Only calculate the qc metrics that are required for the requested filtering.
        // These are computed once, per-sample, on the raw counts:
        //   - total_counts (.obs)      -> used to filter cells on --min_counts/--max_counts
        //   - num_nonzero_vars (.obs)  -> used to filter cells on --min_genes_per_cell/--max_genes_per_cell
        //   - num_nonzero_obs (.var)   -> used to filter genes on --min_cells_per_gene
        // Mitochondrial gene detection is handled by the blocks below.
        def filter_on_counts = state.min_counts != null || state.max_counts != null
        def filter_on_genes_per_cell = state.min_genes_per_cell != null || state.max_genes_per_cell != null
        def filter_on_cells_per_gene = state.min_cells_per_gene != null
        def args = [
          "id": id,
          "input": state.input,
          "top_n_vars": [],
          "output_obs_num_nonzero_vars": filter_on_genes_per_cell ? "num_nonzero_vars" : null,
          "output_obs_total_counts_vars": filter_on_counts ? "total_counts" : null,
          "output_var_num_nonzero_obs": filter_on_cells_per_gene ? "num_nonzero_obs" : null,
          "output_var_total_counts_obs": null,
          "output_var_obs_mean": null,
          "output_var_pct_dropout": null,
          "output": state.output,
          "modality": "rna",
          "layer": state.layer,
        ]
        
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
    // cell filtering on total counts per cell (.obs total_counts, computed by qc above)
    | delimit_counts.run(
      key: "rna_filter_total_counts",
      runIf: { id, state -> state.min_counts != null || state.max_counts != null },
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "obs_count_column": ["total_counts"],
          "obs_name_filter": ["filter_total_counts"],
          "min_count": state.min_counts,
          "max_count": state.max_counts,
        ]
      },
      toState: ["input": "output"]
    )
    // cell filtering on the number of genes per cell (.obs num_nonzero_vars)
    | delimit_counts.run(
      key: "rna_filter_genes_per_cell",
      runIf: { id, state -> state.min_genes_per_cell != null || state.max_genes_per_cell != null },
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "obs_count_column": ["num_nonzero_vars"],
          "obs_name_filter": ["filter_genes_per_cell"],
          "min_count": state.min_genes_per_cell,
          "max_count": state.max_genes_per_cell,
        ]
      },
      toState: ["input": "output"]
    )
    // gene filtering on the number of cells per gene (.var num_nonzero_obs)
    | delimit_counts.run(
      key: "rna_filter_cells_per_gene",
      runIf: { id, state -> state.min_cells_per_gene != null },
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "var_count_column": ["num_nonzero_obs"],
          "var_name_filter": ["filter_cells_per_gene"],
          "min_count": state.min_cells_per_gene,
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      key: "rna_do_filter",
      fromState: {id, state ->
        // do_filter does not need a layer argument because it filters all layers
        // from a modality. Only the filter columns that were actually created
        // (i.e. for which a threshold was requested) are referenced here.
        def obs_filter = []
        if (state.min_counts != null || state.max_counts != null) {
          obs_filter += ["filter_total_counts"]
        }
        if (state.min_genes_per_cell != null || state.max_genes_per_cell != null) {
          obs_filter += ["filter_genes_per_cell"]
        }
        if (state.var_name_mitochondrial_genes) {
          obs_filter += ["filter_mitochondrial"]
        }
        if (state.var_name_ribosomal_genes) {
          obs_filter += ["filter_ribosomal"]
        }
        def var_filter = []
        if (state.min_cells_per_gene != null) {
          var_filter += ["filter_cells_per_gene"]
        }
        def stateMapping = [
          input: state.input,
          // If scrublet is skipped, the output should be set to the workflow output
          output: state.workflow_output
        ]
        // Only pass the filter arguments when columns were actually created;
        // an empty list would be passed as a single empty string to do_filter.
        if (obs_filter) {
          stateMapping["obs_filter"] = obs_filter
        }
        if (var_filter) {
          stateMapping["var_filter"] = var_filter
        }
        return stateMapping
      },
      toState: ["input": "output"]
    )
    // doublet calling
    | filter_with_scrublet.run(
      runIf: { id, state ->
        !state.skip_scrublet_doublet_detection
      },
      fromState: [
        "input": "input",
        "layer": "layer",
        "scrublet_score_threshold": "scrublet_score_threshold",
        "output": "workflow_output"
      ],
      args: [output_compression: "gzip"],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}