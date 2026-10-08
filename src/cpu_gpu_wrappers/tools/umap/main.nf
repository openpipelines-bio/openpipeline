workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | umap_gpu.run(
      runIf: { id, state -> state.device_type =="gpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "uns_neighbors": "uns_neighbors",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "min_dist": "min_dist",
        "spread": "spread",
        "num_components": "num_components",
        "max_iter": "max_iter",
        "alpha": "alpha",
        "negative_sample_rate": "negative_sample_rate",
        "init_pos": "init_pos",
        "random_state": "random_state"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | umap_cpu.run(
      runIf: { id, state -> state.device_type =="cpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "uns_neighbors": "uns_neighbors",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "min_dist": "min_dist",
        "spread": "spread",
        "num_components": "num_components",
        "max_iter": "max_iter",
        "alpha": "alpha",
        "negative_sample_rate": "negative_sample_rate",
        "init_pos": "init_pos",
        "gamma": "gamma"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
