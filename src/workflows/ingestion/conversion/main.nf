workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | runEach(
      components: [from_10xh5_to_h5mu, from_h5ad_to_h5mu, from_10xmtx_to_h5mu],
      filter: { id, state, component ->
        def componentNameMapper = [
          "10xh5": "from_10xh5_to_h5mu",
          "10xmtx": "from_10xmtx_to_h5mu",
          "h5ad": "from_h5ad_to_h5mu"
        ]
        componentNameMapper[state.input_type] == component.config.name
      },
      fromState: { id, state, component ->
        def passed_state = [
          input: state.input,
          compression: "gzip",
          output: state.output
        ]
        if (component.name == "from_h5ad_to_h5mu") {
          passed_state.modality = state.modality
        }
        passed_state
      },
      toState: ["output": "output"]
    )
    | setState(["output": "output"])
  
  emit:
  output_ch
}