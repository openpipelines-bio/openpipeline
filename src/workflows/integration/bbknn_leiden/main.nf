workflow run_wf {
  take:
  input_ch

  main:
  bbknn_ch = input_ch
    // Make sure there is not conflict between the output from this workflow
    // And the output from any of the components
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // compute bbknn graph
    | bbknn.run(
      fromState: { id, state ->
        [
          input: state.input,
          modality: state.modality,
          obsm_input: state.obsm_input,
          obs_batch: state.obs_batch,
          uns_output: state.uns_output,
          obsp_distances: state.obsp_distances,
          obsp_connectivities: state.obsp_connectivities,
          n_neighbors_within_batch: state.n_neighbors_within_batch,
          n_pcs: state.n_pcs,
          n_trim: state.n_trim
        ]
      },
      // use map when viash 0.7.6 is released
      // related to https://github.com/viash-io/viash/pull/515
      // fromState: [
      //   "input": "input",
      //   "obsm_input": "obsm_input",
      //   "obs_batch": "obs_batch",
      //   "modality": "modality",
      //   "uns_output": "uns_output",
      //   "obsp_distances": "obsp_distances",
      //   "obsp_connectivities": "obsp_connectivities",
      //   "n_neighbors_within_batch": "n_neighbors_within_batch",
      //   "n_pcs": "n_pcs",
      //   "n_trim": "n_trim"
      // ],
      toState: [
        "input": "output"
      ]
    )
    | neighbors_leiden_umap.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_integrated",
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
  with_leiden_ch = bbknn_ch
    | filter{id, state -> state.leiden_resolution}
    // run leiden on the bbknn graph
    | leiden.run(
      fromState: [
        "input": "input",
        "obsp_connectivities": "obsp_connectivities",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution",
        "modality": "modality"
      ],
      toState: [
        "input": "output"
      ]
    )
    // move obsm leiden cluster dataframe to obs
    | move_obsm_to_obs.run(
      fromState:
        [
          "input": "input",
          "obsm_key": "obs_cluster",
          "modality": "modality",
        ],
      toState: ["input": "output"]
    )

  without_leiden_ch = bbknn_ch
    | filter{id, state -> !state.leiden_resolution}
  
  output_ch = with_leiden_ch.mix(without_leiden_ch)
    // run umap on the bbknn graph
    | umap.run(
      fromState: { id, state ->
       [
          "input": state.input,
          "uns_neighbors": state.uns_output,
          "obsm_output": state.obsm_umap,
          "modality": state.modality,
          "output": state.workflow_output,
          "output_compression": "gzip"
       ]
      },
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}