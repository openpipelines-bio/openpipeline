methods = [
  densmap,
  diffusion_map,
  ivis,
  lmds,
  neuralee,
  pca,
  phate,
  pymde,
  simlr,
  tsne,
  umap
]

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | openproblems_dr_h5mu_to_h5ad.run(
      fromState: [
        "input",
        "input_modality",
        "input_layer_counts",
        "input_layer_normalized",
        "input_var_hvg_score"
      ],
      toState: [
        "method_input": "output"
      ]
    )

    | runEach(
      components: methods,
      runIf: { id, state, comp ->
        state.method_id == comp.config.name
      },
      fromState: [
        "input": "method_input"
      ],
      toState: [
        "method_output": "output"
      ]
    )

    | openproblems_dr_h5ad_to_h5mu.run(
      fromState: [
        "input_dataset": "input",
        "input_output": "method_output",
        "input_modality": "input_modality",
        "output_obsm_key": "output_obsm_key",
      ],
      toState: [
        "output": "output"
      ]
    )

  emit:
  output_ch
}
