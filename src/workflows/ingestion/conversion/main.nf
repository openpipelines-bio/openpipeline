workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // translate the input_type argument to the component name that needs to be run
    | map { id, state ->
      def componentNameMapper = [
        "10xh5": from_10xh5_to_h5mu,
        "10xmtx": from_10xmtx_to_h5mu,
        "h5ad": from_h5ad_to_h5mu
      ]
      def component = componentNameMapper[state.input_type]
      def new_state = state + ["component": component]
      [id, new_state]
    }
    | runEach(
      components: [from_10xh5_to_h5mu, from_h5ad_to_h5mu, from_10xmtx_to_h5mu],
      filter: { id, state, component ->
        state.component == component
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
      toState: {id, output, state, comp ->
        ["output": output.output]
      },
      auto: [publish: true],
    )
  
  emit:
  output_ch
}