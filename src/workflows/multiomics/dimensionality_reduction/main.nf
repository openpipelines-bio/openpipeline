workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | pca.run(
      fromState: [
        "input": "input", 
        "obsm_output": "obsm_pca",
        "var_input": "var_pca_feature_selection",
        "modality": "modality",
        "overwrite": "pca_overwrite",
        "layer": "layer",
        "varm_output": "pca_loadings_varm_output",
        "uns_output": "pca_variance_uns_output",
      ],
      toState: ["input": "output"]
    )
    | neighbors_leiden_umap.run(
      fromState: [
        "input": "input",
        "obsm_input": "obsm_pca",
        "modality": "modality",
        "uns_neighbors": "uns_neighbors",
        "obsp_neighbor_distances": "obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "obsp_neighbor_connectivities",
        "output": "workflow_output",
        "obsm_umap": "obsm_umap",
      ],
      toState: ["output": "output"],
      args: [
        "leiden_resolution": [] // disable leiden
      ]
    )
    | setState(["output"])

  emit:
  output_ch
}