workflow run_wf {
  take:
  input_ch

  main:
  neighbors_ch = input_ch
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
    | find_neighbors.run(
      fromState: [
        "input": "input",
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_output", // use output from scvi as input for neighbors,
        "modality": "modality"
      ],
      toState: ["input": "output"]
    )

  with_leiden_ch = neighbors_ch
    | filter{id, state -> state.leiden_resolution}
    | leiden.run(
      fromState: [
        "input": "input",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution",
        "modality": "modality",
      ],
      toState: ["input": "output"]
    )
    | move_obsm_to_obs.run(
      fromState: [
          "input": "input",
          "obsm_key": "obs_cluster",
          "modality": "modality",
      ],
      toState: ["input": "output"]
    )


  without_leiden_ch = neighbors_ch
    | filter{id, state -> !state.leiden_resolution}

  output_ch = with_leiden_ch.mix(without_leiden_ch)
    | umap.run(
      fromState: {id, state -> [
        "input": state.input,
        "uns_neighbors": state.uns_neighbors,
        "obsm_output": state.obsm_umap,
        "modality": state.modality,
        "output": state.workflow_output,
        "output_compression": "gzip"
        ]
      },
      auto: [ publish: true ],
      toState: { id, output, state ->
        [ 
          output: output.output, 
          output_model: state.output_model
        ]
      }
    )

  emit:
  output_ch
}
