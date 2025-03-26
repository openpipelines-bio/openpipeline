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
    | setState(["output"])

  emit:
  output_ch
}
