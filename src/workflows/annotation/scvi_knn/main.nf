workflow run_wf {
  take:
    input_ch

  main:
    modalities_ch = input_ch
      // Set aside the output for this workflow to avoid conflicts
      | map {id, state -> 
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
      }
      // Allign query and reference datasets
      | align_query_reference.run(
        fromState: [
          "input": "input",
          "modality": "modality",
          "input_layer": "input_layer",
          "input_layer_lognormalized": "input_layer_lognormalized",
          "input_obs_batch": "input_obs_batch_label",
          "input_var_gene_names": "input_var_gene_names",
          "reference": "reference",
          "reference_layer": "reference_layer",
          "reference_layer_lognormalized": "reference_layer_lognormalized",
          "reference_obs_batch": "reference_obs_batch_label",
          "reference_obs_label": "reference_obs_target",
          "reference_var_gene_names": "reference_var_gene_names",
          "input_reference_gene_overlap": "input_reference_gene_overlap",
          "overwrite_existing_key": "overwrite_existing_key"
        ],
        args: [
          "input_id": "query",
          "reference_id": "reference",
          "output_layer": "_counts",
          "output_layer_lognormalized": "_log_normalized",
          "output_var_gene_names": "_gene_names",
          "output_obs_batch": "_sample_id",
          "output_obs_label": "_cell_type",
          "output_obs_id": "_dataset",
          "output_var_common_genes": "_common_vars",
          "align_layers_lognormalized_counts": "true"
        ],
        toState: [
          "input": "output_query",
          "reference": "output_reference"
        ]
      )
      
      | split_modalities.run(
        fromState: {id, state ->
          def newState = ["input": state.input, "id": id]
        },
        toState: ["output": "output", "output_types": "output_types"]
      )
      | flatMap {id, state ->
        def outputDir = state.output
        def types = readCsv(state.output_types.toUriString())
        
        types.collect{ dat ->
          // def new_id = id + "_" + dat.name
          def new_id = id // it's okay because the channel will get split up anyways
          def new_data = outputDir.resolve(dat.filename)
          [ new_id, state + ["input": new_data, modality: dat.name]]
        }
      }
      // Remove arguments from split modalities from state
      | map {id, state -> 
        def keysToRemove = ["output_types", "output_files"]
        def newState = state.findAll{it.key !in keysToRemove}
        [id, newState]
      }
      | view {"After splitting modalities: $it"}

    rna_ch = modalities_ch

      | filter{id, state -> state.modality == "rna"}

      // Concatenate query and reference datasets prior to integration
      // Only concatenate rna modality in this channel
      | concatenate_h5mu.run(
        fromState: { id, state -> 
          ["input": [state.input, state.reference]]
        },
        args: [
          "input_id": ["query", "reference"],
          "modality": "rna",
          "other_axis_mode": "move"
        ],
        toState: ["input": "output"]
      )

      // Calculate HVG across query and reference
      | highly_variable_features_scanpy.run(
        fromState: [
          "input": "input",
          "modality": "modality",
          "n_top_features": "n_hvg"
        ],
        args: [
          "layer": "_log_normalized",
          "var_input": "_common_vars",
          "var_name_filter": "_common_hvg",
          "obs_batch_key": "_sample_id"
        ],
        toState: [
          "input": "output"
        ]
      )
      | delete_layer.run(
        key: "delete_aligned_lognormalized_counts_layer",
        fromState: [
          "input": "input",
          "modality": "modality",
        ],
        args: [
          "layer": "_log_normalized",
          "missing_ok": "true"
        ],
        toState: [
          "input": "output"
        ]
      )
      // Run scvi integration with leiden clustering
      | scvi_leiden_workflow.run(
        fromState: { id, state -> [
          "id": id,
          "input": state.input,
          "modality": state.modality,
          "obsm_output": state.output_obsm_integrated,
          "leiden_resolution": state.leiden_resolution,
          "early_stopping": state.scvi_early_stopping,
          "early_stopping_monitor": state.scvi_early_stopping_monitor,
          "early_stoping_patience": state.scvi_early_stoping_patience,
          "early_stopping_min_delta": state.scvi_early_stopping_min_delta,
          "max_epochs": state.scvi_max_epochs,
          "reduce_lr_on_plateau": state.scvi_reduce_lr_on_plateau,
          "lr_factor": state.scvi_lr_factor,
          "lr_patience": state.scvi_lr_patience,
          "sanitize_ensembl_ids": state.sanitize_ensembl_ids
        ]},
        args: [
          "var_input": "_common_hvg",
          "uns_neighbors": "scvi_integration_neighbors",
          "obsp_neighbor_distances": "scvi_integration_distances",
          "obsp_neighbor_connectivities": "scvi_integration_connectivities",
          "obs_cluster": "scvi_integration_leiden",
          "obsm_umap": "X_leiden_scvi_umap",
          "obs_batch": "_sample_id",
          "layer": "_counts"
        ],
        toState: ["input": "output"]
      )
      | delete_layer.run(
        key: "delete_aligned_counts_layer",
        fromState: [
          "input": "input",
          "modality": "modality",
        ],
        args: [
          "layer": "_counts",
          "missing_ok": "true"
        ],
        toState: [
          "input": "output"
        ]
      )
      // Split integrated dataset back into a separate reference and query dataset
      | split_h5mu.run(
        fromState: [
          "input": "input",
          "modality": "modality"
        ],
        args: [
          "obs_feature": "_dataset",
          "output_files": "sample_files.csv",
          "drop_obs_nan": "true",
          "output": "ref_query"
        ],
        toState: [ 
          "output": "output", 
          "output_files": "output_files" 
        ]
      )
      // map the integrated query and reference datasets back to the state
      | map {id, state ->
          def outputDir = state.output
          if (workflow.stubRun) {
            def output_files = outputDir.listFiles()
            def new_state = state + [
                "input": output_files[0],
                "reference": output_files[1],
            ]
            return [id, new_state]
          }
          def files = readCsv(state.output_files.toUriString())
          def query_file = files.findAll{ dat -> dat.name == 'query' }
          assert query_file.size() == 1, 'there should only be one query file'
          def reference_file = files.findAll{ dat -> dat.name == 'reference' }
          assert reference_file.size() == 1, 'there should only be one reference file'
          def integrated_query = outputDir.resolve(query_file.filename)
          def integrated_reference = outputDir.resolve(reference_file.filename)
          def newKeys = ["input": integrated_query, "reference": integrated_reference]
          [id, state + newKeys]
      }
      // remove keys from split files
      | map {id, state -> 
        def keysToRemove = ["output_files"]
        def newState = state.findAll{it.key !in keysToRemove}
        [id, newState]
      }

      // Perform KNN label transfer from reference to query
      | knn.run(
        fromState: [
           "input": "input",
           "modality": "modality",
           "input_obsm_features": "output_obsm_integrated",
           "reference": "reference",
           "reference_obsm_features": "output_obsm_integrated",
           "reference_obs_targets": "reference_obs_target",
           "output_obs_predictions": "output_obs_predictions",
           "output_obs_probability": "output_obs_probability",
           "weights": "knn_weights",
           "n_neighbors": "knn_n_neighbors",
           "output": "workflow_output"
        ],
        args:[
          "output_compression": "gzip"
        ],
        toState: ["input": "output"]
        // toState: {id, output, state -> ["output": output.output]}
      )
      | view {"After processing RNA modality: $it"}


    other_mod_ch = modalities_ch
      | filter{id, state -> state.modality != "rna"}

    output_ch = rna_ch.mix(other_mod_ch)
      | groupTuple(by: 0, sort: "hash")
      | map { id, states ->
          def new_input = states.collect{it.input}
          def modalities = states.collect{it.modality}.unique()
          def other_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
          }.minus(["output", "input", "modality", "reference"])
          def new_state = other_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "Arguments should be the same across modalities. Please report this \
                                                 as a bug. Argument name: $argument_name, \
                                                 argument value: $argument_values"
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
          }
          [id, new_state + ["input": new_input, "modalities": modalities]]
      }
      | merge.run(
        fromState: ["input": "input"],
        toState: ["output": "output"],
      )
      | setState(["output"])
    
  emit:
    output_ch
}
