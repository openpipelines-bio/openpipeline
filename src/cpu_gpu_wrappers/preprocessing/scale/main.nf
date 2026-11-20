workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | scale_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "input_layer": "input_layer",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "zero_center": "zero_center",
        "max_value": "max_value"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | scale_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "input_layer": "input_layer",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "zero_center": "zero_center",
        "max_value": "max_value"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
