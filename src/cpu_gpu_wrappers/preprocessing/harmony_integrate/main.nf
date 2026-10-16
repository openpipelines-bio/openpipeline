workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // -- GPU variant (rapids-singlecell) --
    | harmony_integrate_gpu.run(
      runIf: { id, state -> state.device_type == "gpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "obs_covariates": "obs_covariates",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "overwrite": "overwrite",
        "theta": "theta",
        "flavor": "flavor",
        "n_clusters": "n_clusters",
        "max_iter_harmony": "max_iter_harmony",
        "random_state": "random_state"
      ],
      toState: ["output": "output"]
    )
    // -- CPU variant (openpipeline) --
    | harmony_integrate_cpu.run(
      runIf: { id, state -> state.device_type == "cpu" },
      fromState: [
        "input": "input",
        "output": "output",
        "modality": "modality",
        "obsm_input": "obsm_input",
        "obs_covariates": "obs_covariates",
        "obsm_output": "obsm_output",
        "output_compression": "output_compression",
        "theta": "theta"
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
