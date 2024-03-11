workflow run_wf {
  take:
  input_ch

  main:
  neighbors_ch = input_ch
    // Make sure there is not conflict between the output from this workflow
    // And the output from any of the components
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // run harmonypy
    | harmonypy.run(
      fromState: [
          "input": "input",
          "modality": "modality",
          "obsm_input": "embedding",
          "obs_covariates": "obs_covariates",
          "obsm_output": "obsm_integrated",
          "theta": "theta"
      ],
      toState: ["input": "output"]
    )
    
    // run knn
    | find_neighbors.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_integrated"
      ],
      toState: ["input": "output"]
    )

  with_leiden_ch = neighbors_ch
    | filter{id, state -> state.leiden_resolution}
    // run leiden clustering
    | leiden.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution"
      ],
      toState: ["input": "output"]
    )
    // move obsm to obs
    | move_obsm_to_obs.run(
      fromState: 
        [
          "input": "input",
          "obsm_key": "obs_cluster",
          "modality": "modality",
        ],
      toState: ["input": "output"]
    )

  without_leiden_ch = neighbors_ch
    | filter{id, state -> !state.leiden_resolution}

  output_ch = with_leiden_ch.mix(without_leiden_ch)
    // run umap
    | umap.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": state.modality,
          "obsm_input": state.obsm_integrated,
          "obsm_output": state.obsm_umap,
          "uns_neighbors": state.uns_neighbors,
          "output": state.workflow_output,
          "output_compression": "gzip"
        ]
      },
      toState: { id, output, state ->
        [ output: output.output ]
      },
      auto: [ publish: true ]
    )

  emit:
  output_ch
}
