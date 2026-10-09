workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | neighbors_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "output_compression": "output_compression",
        "uns_output": "uns_output",
        "obsp_distances": "obsp_distances",
        "obsp_connectivities": "obsp_connectivities",
        "num_neighbors": "num_neighbors",
        "metric": "metric",
        "n_pcs": "n_pcs",
        "algorithm": "algorithm",
        "method": "method",
        "random_state": "random_state"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | neighbors_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "output_compression": "output_compression",
        "uns_output": "uns_output",
        "obsp_distances": "obsp_distances",
        "obsp_connectivities": "obsp_connectivities",
        "num_neighbors": "num_neighbors",
        "metric": "metric",
        "seed": "random_state"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
