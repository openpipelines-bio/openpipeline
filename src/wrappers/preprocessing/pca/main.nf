workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | pca_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "layer": "layer",
        "var_input": "var_input",
        "num_components": "num_components",
        "chunked": "chunked",
        "chunk_size": "chunk_size",
        "random_state": "random_state",
        "output_compression": "output_compression",
        "obsm_output": "obsm_output",
        "varm_output": "varm_output",
        "uns_output": "uns_output",
        "overwrite": "overwrite"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | pca_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "layer": "layer",
        "var_input": "var_input",
        "num_components": "num_components",
        "chunked": "chunked",
        "chunk_size": "chunk_size",
        "seed": "random_state",
        "output_compression": "output_compression",
        "obsm_output": "obsm_output",
        "varm_output": "varm_output",
        "uns_output": "uns_output",
        "overwrite": "overwrite"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
