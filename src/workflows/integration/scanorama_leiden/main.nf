workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
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
    | setState(["output"])

  emit:
  output_ch
}
