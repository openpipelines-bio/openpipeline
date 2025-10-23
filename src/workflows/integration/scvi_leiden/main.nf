workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | map {id, state -> 
      assert state.input || state.tiledb_input_uri: "Either --input or --tiledb_input_uri must be defined"
      assert !(state.input && state.tiledb_input_uri): "Values were specified for both --input and --tiledb_input_uri. Please choose one."
      assert state.tiledb_s3_region || !state.tiledb_input_uri: "Specifying 'tiledb_s3_region' also requires 'tiledb_input_uri'."
      assert !state.tiledb_input_uri || state.layer: "When using tileDB input, you must specify a layer using --layer"
      [id, state]
    }
    | from_tiledb_to_h5mu.run(
      runIf: {id, state -> state.tiledb_input_uri},
      fromState: { id, state -> 
        [
          "input_uri": state.tiledb_input_uri,
          "s3_region": state.tiledb_s3_region,
          "endpoint": state.tiledb_endpoint,
          "input_modality": state.modality,
          "output_modality": state.modality,
          "input_layers": [ state.layer ],
          "s3_no_sign_request": state.tiledb_s3_no_sign_request
        ]
      },
      toState: [
        "input": "output",
      ]
    )
    | scvi.run(
      fromState: [
        "input": "input",
        "obs_batch": "obs_batch",
        "obs_size_factor": "obs_size_factor",
        "obs_categorical_covariate": "obs_categorical_covariate",
        "obs_continuous_covariate": "obs_continuous_covariate",
        "obsm_output": "obsm_output",
        "var_input": "var_input",
        "early_stopping": "early_stopping",
        "early_stopping_monitor": "early_stopping_monitor",
        "early_stopping_patience": "early_stopping_patience",
        "early_stopping_min_delta": "early_stopping_min_delta",
        "max_epochs": "max_epochs",
        "reduce_lr_on_plateau": "reduce_lr_on_plateau",
        "lr_factor": "lr_factor",
        "lr_patience": "lr_patience",
        "output_model": "output_model",
        "modality": "modality",
        "input_layer": "layer",
      ],
      toState: [
        "input": "output", 
        "output_model": "output_model"
      ],
    )
    | neighbors_leiden_umap.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_output",
        "output": "workflow_output",
        "uns_neighbors": "uns_neighbors",
        "obsp_neighbor_distances": "obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "obsp_neighbor_connectivities",
        "leiden_resolution": "leiden_resolution",
        "obs_cluster": "obs_cluster",
        "obsm_umap": "obsm_umap",
      ],
      toState: ["output": "output"]
    )
    | move_mudata_obsm_to_tiledb.run(
      runIf: {id, state -> state.tiledb_input_uri},
      fromState: {id, state ->
        [
          "input_mudata": state.output,
          "modality": state.modality,
          "input_uri": state.tiledb_input_uri,
          "s3_region": state.tiledb_s3_region,
          "endpoint": state.tiledb_endpoint,
          "output_modality": state.modality,
          "obsm_input": [state.obsm_umap, state.obsm_output],
          "output_tiledb": state.output_tiledb,
          "s3_no_sign_request": state.tiledb_s3_no_sign_request
        ]
      },
      toState: ["previous_output_tiledb": "output_tiledb"]
    )
    | move_mudata_obs_to_tiledb.run(
      runIf: {id, state -> state.tiledb_input_uri},
      fromState: {id, state ->
        def new_state = [
          "input_mudata": state.output,
          "modality": state.modality,
          "s3_region": state.tiledb_s3_region,
          "endpoint": state.tiledb_endpoint,
          "output_modality": state.modality,
          "obs_input": state.leiden_resolution.collect{state.obs_cluster + "_" + it},
          "output_tiledb": state.output_tiledb,
          "s3_no_sign_request": state.tiledb_s3_no_sign_request
        ]
        if (state.previous_output_tiledb && file(state.previous_output_tiledb).exists()) {
          new_state += ["input_dir": state.previous_output_tiledb]
        } else {
          new_state += ["input_uri": state.tiledb_input_uri]
        }
        return new_state
      },
      toState: ["output_tiledb": "output_tiledb"]
    )
    | setState(["output", "output_model", "output_tiledb"])

  emit:
  output_ch
}
