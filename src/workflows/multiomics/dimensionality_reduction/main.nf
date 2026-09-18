workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // The GPU (rapids-singlecell) PCA errors out on genes with zero expression,
    // which upstream cell filtering can leave behind. The CPU PCA tolerates them,
    // so only the GPU path needs this.
    | filter_genes.run(
      runIf: {id, state -> state.device_type == "gpu"},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": state.modality,
          "layer": state.layer,
        ]
      },
      args: ["min_counts": 1],
      toState: ["input": "output"]
    )
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
        "device_type": "device_type",
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
        "device_type": "device_type",
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