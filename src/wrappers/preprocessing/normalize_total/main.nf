workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | normalize_total_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "input_layer": "input_layer",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "target_sum": "target_sum",
        "exclude_highly_expressed": "exclude_highly_expressed",
        "max_fraction": "max_fraction"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | normalize_total_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "input_layer": "input_layer",
        "output_layer": "output_layer",
        "output_compression": "output_compression",
        "target_sum": "target_sum",
        "exclude_highly_expressed": "exclude_highly_expressed"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
