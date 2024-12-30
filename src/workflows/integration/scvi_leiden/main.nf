workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | scvi.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "obs_batch": state.obs_batch,
          "obsm_output": state.obsm_output,
          "var_input": state.var_input,
          "early_stopping": state.early_stopping,
          "early_stopping_monitor": state.early_stopping_monitor,
          "early_stopping_patience": state.early_stopping_patience,
          "early_stopping_min_delta": state.early_stopping_min_delta,
          "max_epochs": state.max_epochs,
          "reduce_lr_on_plateau": state.reduce_lr_on_plateau,
          "lr_factor": state.lr_factor,
          "lr_patience": state.lr_patience,
          "output_model": state.output_model,
          "modality": state.modality,
          "input_layer": state.layer,
       ]
      },
    // use map when viash 0.7.6 is released
    // related to https://github.com/viash-io/viash/pull/515
    //   fromState: [
    //     "input": "input",
    //     "obs_batch": "obs_batch",
    //     "obsm_output": "obsm_output",
    //     "var_input": "var_input",
    //     "early_stopping": "early_stopping",
    //     "early_stopping_monitor": "early_stopping_monitor",
    //     "early_stopping_patience": "early_stopping_patience",
    //     "early_stopping_min_delta": "early_stopping_min_delta",
    //     "max_epochs": "max_epochs",
    //     "reduce_lr_on_plateau": "reduce_lr_on_plateau",
    //     "lr_factor": "lr_factor",
    //     "lr_patience": "lr_patience",
    //     "output_model": "output_model",
    //     "modality": "modality"
    //   ],
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
    | setState(["output", "output_model"])

  emit:
  output_ch
}
