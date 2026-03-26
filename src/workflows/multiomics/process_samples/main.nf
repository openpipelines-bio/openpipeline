workflow run_wf {
  take:
    input_ch

  main:
    singlesample_ch = input_ch
      // Make sure there is not conflict between the output from this workflow
      // And the output from any of the components
      | map {id, state ->
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
      }

      // If requested to be detected, make sure the mitochondrial and ribosomal genes
      // are added to the input of the qc metrics calculation
      | map {id, state ->
        def var_qc_default = [state.highly_variable_features_var_output]
        if (state.var_name_mitochondrial_genes) {
          var_qc_default.add(state.var_name_mitochondrial_genes)
        }
        if (state.var_name_ribosomal_genes) {
          var_qc_default.add(state.var_name_ribosomal_genes)
        }
        def newState = state + ["var_qc_metrics": var_qc_default.join(",")]
        [id, newState]
      }

      | process_singlesample.run(
        fromState: {id, state ->
          [
            "id": id,
            "input": state.input,
            "rna_layer": state.rna_layer,
            "prot_layer": state.prot_layer,
            "gdo_layer": state.gdo_layer,
            "add_id_to_obs": state.add_id_to_obs,
            "add_id_obs_output": state.add_id_obs_output,
            "add_id_make_observation_keys_unique": state.add_id_make_observation_keys_unique,
            "rna_min_counts": state.rna_min_counts,
            "rna_max_counts": state.rna_max_counts,
            "rna_min_percentile_counts": state.rna_min_percentile_counts,
            "rna_max_percentile_counts": state.rna_max_percentile_counts,
            "rna_log_transform_total_counts": state.rna_log_transform_total_counts,
            "rna_min_genes_per_cell": state.rna_min_genes_per_cell,
            "rna_max_genes_per_cell": state.rna_max_genes_per_cell,
            "rna_min_cells_per_gene": state.rna_min_cells_per_gene,
            "rna_min_fraction_mito": state.rna_min_fraction_mito,
            "rna_max_fraction_mito": state.rna_max_fraction_mito,
            "rna_min_fraction_ribo": state.rna_min_fraction_ribo,
            "rna_max_fraction_ribo": state.rna_max_fraction_ribo,
            "skip_scrublet_doublet_detection": state.skip_scrublet_doublet_detection,
            "prot_min_counts": state.prot_min_counts,
            "prot_max_counts": state.prot_max_counts,
            "prot_min_percentile_counts": state.prot_min_percentile_counts,
            "prot_max_percentile_counts": state.prot_max_percentile_counts,
            "prot_log_transform_total_counts": state.prot_log_transform_total_counts,
            "prot_min_proteins_per_cell": state.prot_min_proteins_per_cell,
            "prot_max_proteins_per_cell": state.prot_max_proteins_per_cell,
            "prot_min_cells_per_protein": state.prot_min_cells_per_protein,
            "gdo_min_counts": state.gdo_min_counts,
            "gdo_max_counts": state.gdo_max_counts,
            "gdo_min_guides_per_cell": state.gdo_min_guides_per_cell,
            "gdo_max_guides_per_cell": state.gdo_max_guides_per_cell,
            "gdo_min_cells_per_guide": state.gdo_min_cells_per_guide,
            "var_gene_names": state.var_gene_names,
            "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
            "obs_name_mitochondrial_fraction": state.obs_name_mitochondrial_fraction,
            "mitochondrial_gene_regex": state.mitochondrial_gene_regex,
            "var_name_ribosomal_genes": state.var_name_ribosomal_genes,
            "obs_name_ribosomal_fraction": state.obs_name_ribosomal_fraction,
            "ribosomal_gene_regex": state.ribosomal_gene_regex
            ]
        },
        toState: ["input": "output"]
      )
  
    def singlesample_arguments = [
      "rna_min_counts",
      "rna_max_counts",
      "rna_min_percentile_counts",
      "rna_max_percentile_counts",
      "rna_min_genes_per_cell",
      "rna_max_genes_per_cell",
      "rna_min_cells_per_gene",
      "rna_min_fraction_mito",
      "rna_max_fraction_mito",
      "rna_min_fraction_ribo",
      "rna_max_fraction_ribo",
      "skip_scrublet_doublet_detection",
      "prot_min_counts",
      "prot_max_counts",
      "prot_min_percentile_counts",
      "prot_max_percentile_counts",
      "prot_min_proteins_per_cell",
      "prot_max_proteins_per_cell",
      "prot_min_cells_per_protein",
      "gdo_min_counts",
      "gdo_max_counts",
      "gdo_min_guides_per_cell",
      "gdo_max_guides_per_cell",
      "gdo_min_cells_per_guide",
      "var_gene_names",
      "var_name_mitochondrial_genes",
      "obs_name_mitochondrial_fraction",
      "mitochondrial_gene_regex",
      "var_name_ribosomal_genes",
      "obs_name_ribosomal_fraction",
      "ribosomal_gene_regex"
    ]

    concat_ch = singlesample_ch
      // Collect all samples before concatenation
      | toList()

      // Concatenation: single, multi-modal MuData objects into a single multi-modal MuData object. 
      | map { sample_list ->
        def old_ids = sample_list.collect { id, state -> id }
        def states = sample_list.collect { id, state -> state }
        def new_id = "merged"

        // keys in the new state that should not have a unique value across samples
        def new_state_non_unique_values = [
          "input": states.collect{it.input},
          "input_id": old_ids,
          "_meta": ["join_id": old_ids[0]]
        ]
        // Gather the keys from the different states,
        // one state might contain more keys compared to another (so create a set)
        def all_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
        }.minus(["output", "input_id", "input", "_meta", "id"] + singlesample_arguments)
        // Create the new state from the keys, values should be the same across samples
        def new_state = all_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1,
            "Arguments should be the same across samples. Argument name: $argument_name, argument value: $argument_values"
            // take the unique value from the set (there is only one)
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
        }
        def final_state = new_state_non_unique_values + new_state
        [new_id, final_state]
      }
      | concatenate_h5mu.run(
        fromState: [
          "input": "input",
          "input_id": "input_id"
        ],
        toState: {id, output, state -> 
          def keysToRemove = ["input_id"]
          def newState = state.findAll{it.key !in keysToRemove}
          newState + ["input": output.output]
        }, 
      )

      | view {"After concatenation: $it"}

    multisample_ch = concat_ch
      | process_batches.run(
        fromState: {id, state ->
          [
            "id": id,
            "input": state.input,
            "output": state.workflow_output,
            "highly_variable_features_var_output": state.highly_variable_features_var_output,
            "highly_variable_features_obs_batch_key": state.highly_variable_features_obs_batch_key,
            "var_qc_metrics": state.var_qc_metrics,
            "top_n_vars": state.top_n_vars, 
            "pca_overwrite": state.pca_overwrite,
            "rna_layer": state.rna_layer,
            "prot_layer": state.prot_layer,
            "clr_axis": state.clr_axis,
            "rna_enable_scaling": state.rna_enable_scaling,
            "rna_scaling_output_layer": state.rna_scaling_output_layer,
            "rna_scaling_pca_obsm_output": state.rna_scaling_pca_obsm_output,
            "rna_scaling_pca_loadings_varm_output": state.rna_scaling_pca_loadings_varm_output,
            "rna_scaling_pca_variance_uns_output": state.rna_scaling_pca_variance_uns_output,
            "rna_scaling_umap_obsm_output": state.rna_scaling_umap_obsm_output,
            "rna_scaling_max_value": state.rna_scaling_max_value,
            "rna_scaling_zero_center": state.rna_scaling_zero_center,
          ]
        },
        toState: {id, output, state -> 
          [
            "output": output.output,
            "_meta": state._meta,
          ]
        }
      )
      | view {"After process_batches: $it"}

  emit:
    multisample_ch
}
