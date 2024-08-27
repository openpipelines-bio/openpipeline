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
    | scanorama.run(
      fromState: [
        "input": "input",
        "obsm_input": "obsm_input",
        "obs_batch": "obs_batch",
        "obsm_output": "obsm_output",
        "modality": "modality",
        "batch_size": "batch_size",
        "sigma": "sigma",
        "approx": "approx",
        "alpha": "alpha",
        "knn": "knn",
      ],
      toState: ["input": "output"]
    )
    | find_neighbors.run(
      fromState: [
        "input": "input",
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_output",
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
        "modality": "modality"
      ],
      toState: ["input": "output"]
    )
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
    | umap.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "uns_neighbors": state.uns_neighbors,
          "obsm_output": state.obsm_umap,
          "modality": state.modality,
          "output": state.output,
          "output_compression": "gzip"
        ]
      },
      auto: [ publish: true ],
      toState: { id, output, state ->
        [ output: output.output ]
      }
    )

  emit:
  output_ch
}
