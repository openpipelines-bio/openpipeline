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
      fromState: [
          "input": "input",
          "obs_batch": "obs_batch",
          "obsm_output": "obsm_output",
          "var_input": "var_input",
          "var_input_gene_names": "var_input_gene_names",
          "unknown_celltype": "unknown_celltype",
          "scvi_reference_model": "scvi_reference_model",
          "early_stopping": "early_stopping",
          "early_stopping_monitor": "early_stopping_monitor",
          "early_stopping_patience": "early_stopping_patience",
          "early_stopping_min_delta": "early_stopping_min_delta",
          "max_epochs": "max_epochs",
          "reduce_lr_on_plateau": "reduce_lr_on_plateau",
          "lr_factor": "lr_factor",
          "lr_patience": "lr_patience",
          "output_model": "output_model",
          "modality": "modality",
          "input_layer": "layer",
       ],
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
      fromState: [
        "input": "input",
        "uns_neighbors": "uns_neighbors",
        "obsm_output": "obsm_umap",
        "modality": "modality",
        "output": "workflow_output",
        ],
      args: ["output_compression": "gzip"],
      toState: ["output": "output"] 
    )
    | setState(["output", "output_model"])

  emit:
  output_ch
}
