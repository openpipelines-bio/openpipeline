// beyond/atlas_building - per-sample h5mu -> multi-donor atlas with obs["subpopulation"]
//
// Channel convention throughout: [ id, state_map ]

workflow run_wf {
  take:
  input_ch  // [ sample_id, { input, ..params.. } ]

  main:

  // -- 1. Per-sample: add participant ID to obs, then QC/filter/scrublet ------
  persample_ch = input_ch
    | map { id, state ->
        [ id, state + ["workflow_output": state.output] ]
      }
    | add_id.run(
        fromState: { id, state -> [
          "input":                       state.input,
          "input_id":                    id,
          "obs_output":                  state.obs_participant_id,
          "make_observation_keys_unique": true,
          "output_compression":          "gzip",
        ]},
        toState: [ "input": "output" ]
      )
    | rna_singlesample.run(
        fromState: { id, state -> [
          "id":                             id,
          "input":                          state.input,
          "layer":                          state.layer,
          "min_counts":                     state.min_counts,
          "max_counts":                     state.max_counts,
          "min_genes_per_cell":             state.min_genes_per_cell,
          "max_genes_per_cell":             state.max_genes_per_cell,
          "min_cells_per_gene":             state.min_cells_per_gene,
          "max_fraction_mito":              state.max_fraction_mito,
          "obs_name_mitochondrial_fraction":state.obs_name_mitochondrial_fraction,
          "var_name_mitochondrial_genes":   state.var_name_mitochondrial_genes,
          "var_gene_names":                 state.var_gene_names,
          "skip_scrublet_doublet_detection":state.skip_scrublet_doublet_detection,
          "output":                         state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 2. Fan-in: concatenate all per-sample h5mu into one object -------------
  concat_ch = persample_ch
    | map { id, state -> [ "all_samples", id, state ] }
    | groupTuple(by: 0, sort: "hash")
    | map { _, ids, states ->
        [
          "beyond_atlas",
          [
            "input":               states.collect { it.input },
            "input_id":            ids,
            "workflow_output":     states[0].workflow_output,
            "_meta":               [ "join_id": ids[0] ],
            // carry forward parameters needed in later steps
            "layer":               states[0].layer,
            "highly_variable_features_var_output":     states[0].highly_variable_features_var_output,
            "highly_variable_features_obs_batch_key":  states[0].highly_variable_features_obs_batch_key,
            "highly_variable_features_n_top_features": states[0].highly_variable_features_n_top_features,
            "obs_covariates":      states[0].obs_covariates,
            "leiden_resolution":   states[0].leiden_resolution,
            "obs_cluster":         states[0].obs_cluster,
            "model":               states[0].model,
            "reference":           states[0].reference,
            "reference_obs_target":states[0].reference_obs_target,
            "input_var_gene_names":states[0].input_var_gene_names,
            "majority_voting":     states[0].majority_voting,
            "obs_celltypist_pred": states[0].obs_celltypist_pred,
            "subpop_leiden_resolution": states[0].subpop_leiden_resolution,
            "obs_subpopulation":   states[0].obs_subpopulation,
            "obs_participant_id":  states[0].obs_participant_id,
          ]
        ]
      }
    | concatenate_h5mu.run(
        key: "concatenate_per_sample",
        fromState: [
          "input":    "input",
          "input_id": "input_id",
        ],
        toState: { id, output, state ->
          state.findAll { it.key != "input_id" } + [ "input": output.output ]
        }
      )

  // -- 3. Multi-sample normalization + HVF ------------------------------------
  multisample_ch = concat_ch
    | rna_multisample.run(
        fromState: { id, state -> [
          "id":                                      id,
          "input":                                   state.input,
          "layer":                                   state.layer,
          "highly_variable_features_var_output":     state.highly_variable_features_var_output,
          "highly_variable_features_obs_batch_key":  state.highly_variable_features_obs_batch_key,
          "highly_variable_features_n_top_features": state.highly_variable_features_n_top_features,
          "output":                                  state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 4. PCA -----------------------------------------------------------------
  pca_ch = multisample_ch
    | pca.run(
        key: "pca_atlas",
        args: [ "obsm_output": "X_pca", "output_compression": "gzip" ],
        fromState: { id, state -> [
          "input":     state.input,
          "var_input": state.highly_variable_features_var_output,
          "output":    state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 5. Harmony integration -> neighbors -> Leiden -> UMAP --------------------
  integrated_ch = pca_ch
    | harmony_leiden.run(
        fromState: { id, state -> [
          "input":           state.input,
          "embedding":       "X_pca",
          "obs_covariates":  state.obs_covariates,
          "leiden_resolution": state.leiden_resolution,
          "obs_cluster":     state.obs_cluster,
          "output":          state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 6. CellTypist annotation -----------------------------------------------
  annotated_ch = integrated_ch
    | celltypist_workflow.run(
        fromState: { id, state -> [
          "input":                      state.input,
          "model":                      state.model,
          "reference":                  state.reference,
          "reference_obs_target":       state.reference_obs_target,
          "input_var_gene_names":       state.input_var_gene_names,
          "majority_voting":            state.majority_voting,
          "output_obs_predictions":     state.obs_celltypist_pred,
          "output":                     state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 7. Split atlas by broad cell type (obs["celltypist_pred"]) -------------
  split_ch = annotated_ch
    | split_h5mu.run(
        args: [ "output": "split_by_celltype", "output_files": "split_files.csv", "output_compression": "gzip" ],
        fromState: { id, state -> [
          "input":       state.input,
          "obs_feature": state.obs_celltypist_pred,
        ]},
        toState: { id, output, state ->
          state + [
            "split_output_dir":   output.output,
            "split_output_files": output.output_files,
          ]
        }
      )
    // Expand one event per cell type
    | flatMap { id, state ->
        def lines  = state.split_output_files.readLines()
        def header = lines[0].split(",")*.trim()
        lines.drop(1).findAll { !it.startsWith("#") && !it.isBlank() }.collect { line ->
          def values = line.split(",")*.trim()
          def entry  = [header, values].transpose().collectEntries()
          def subtype_id = id + "_" + entry.name
          [ subtype_id,
            state + [
              "input":     state.split_output_dir.resolve(entry.filename),
              "cell_type": entry.name,
              "_meta":     [ "join_id": state._meta?.join_id ?: id ],
            ]
          ]
        }
      }

  // -- 8. Per-cell-type: PCA -> neighbors -> Leiden (subpopulations) ------------
  pertype_ch = split_ch
    | pca.run(
        key: "pca_celltype",
        args: [ "obsm_output": "X_pca_celltype", "output_compression": "gzip", "overwrite": true ],
        fromState: { id, state -> [
          "input":  state.input,
          "output": state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )
    | neighbors_leiden_umap.run(
        fromState: { id, state -> [
          "input":                      state.input,
          "obsm_input":                 "X_pca_celltype",
          "leiden_resolution":          state.subpop_leiden_resolution,
          "obs_cluster":                state.obs_subpopulation,
          "output":                     state.workflow_output,
        ]},
        args: [
          "uns_neighbors":              "neighbors",
          "obsp_neighbor_distances":    "distances",
          "obsp_neighbor_connectivities": "connectivities",
        ],
        toState: [ "input": "output" ]
      )

  // -- 9. Fan-in: concatenate all cell-type objects -> final atlas -------------
  output_ch = pertype_ch
    | map { id, state -> [ state._meta.join_id, id, state ] }
    | groupTuple(by: 0, sort: "hash")
    | map { atlas_id, cell_type_ids, states ->
        [
          atlas_id,
          [
            "input":    states.collect { it.input },
            "input_id": cell_type_ids,
            "output":   states[0].workflow_output,
          ]
        ]
      }
    | concatenate_h5mu.run(
        key: "concatenate_final",
        fromState: [
          "input":    "input",
          "input_id": "input_id",
          "output":   "output",
        ],
        toState: [ "output": "output" ]
      )
    | setState(["output"])

  emit:
  output_ch
}
