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
    | find_neighbors.run(
      fromState: [
        "input": "input",
        "obsm_input": "obsm_pca",
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "modality": "modality",
        "layer": "layer",
      ],
      toState: ["input": "output"]
    )
    | umap.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "uns_neighbors": state.uns_neighbors,
          "output": state.workflow_output,
          "obsm_output": state.obsm_umap,
          "modality": state.modality,
          "output_compression": "gzip"
        ]
      },
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}