workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | leiden_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "obsp_connectivities": "obsp_connectivities",
        "obsm_name": "obsm_name",
        "output_compression": "output_compression",
        "resolution": "resolution",
        "n_iterations": "n_iterations",
        "random_state": "random_state",
        "theta": "theta",
        "use_weights": "use_weights"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | leiden_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "obsp_connectivities": "obsp_connectivities",
        "obsm_name": "obsm_name",
        "output_compression": "output_compression",
        "resolution": "resolution",
        "n_iterations": "n_iterations",
        "seed": "random_state",
        "flavor": "flavor"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
