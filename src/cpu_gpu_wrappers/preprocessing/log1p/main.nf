workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | log1p_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "input_layer": "input_layer",
        "input_obsm": "input_obsm",
        "output_layer": "output_layer",
        "output_obsm": "output_obsm",
        "output_compression": "output_compression",
        "base": "base"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | log1p_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "input_layer": "input_layer",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "base": "base"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
