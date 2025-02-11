workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state ->
      assert (state.leiden_resolution.isEmpty() || state.obs_cluster?.trim()):
        "When leiden_resolution is set, obs_cluster must also be defined."
      [id, state]
    }
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | find_neighbors.run(
      fromState: [
        "input": "input",
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_input",
        "output": "workflow_output",
        "modality": "modality"
      ],
      toState: ["input": "output"]
    )
    | leiden.run(
      runIf: {id, state -> state.leiden_resolution},
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
      runIf: {id, state -> state.leiden_resolution},
      fromState: [
        "input": "input",
        "output": "workflow_output",
        "obsm_key": "obs_cluster",
        "modality": "modality",
      ],
      args: ["output_compression": "gzip"],
      toState: ["input": "output"]
    )
    | umap.run(
      runIf: {id, state -> !state.obsm_umap?.trim()?.isEmpty()},
      fromState: [
          "input": "input",
          "output": "workflow_output",
          "uns_neighbors": "uns_neighbors",
          "obsm_output": "obsm_umap",
          "modality": "modality",
        ],
      args: ["output_compression": "gzip"],
      toState: ["input": "output"]
    )
    | setState(["output": "input"])

  emit:
  output_ch
}