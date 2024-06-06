workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Make sure there is not conflict between the output from this workflow
    // And the output from any of the components
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // download the reference h5ad file
    | download_file.run(
      fromState: { id, state ->
        [
          "input": state.reference_url,
          "output": "reference.h5ad",
          "verbose": "true",
        ]
      },
      toState: [
        "input": "output",
      ]
    )
    // convert the reference h5ad file to h5mu
    | from_h5ad_to_h5mu.run(
        fromState: { id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "output": "reference.h5mu",
        ]
      },
      toState: { id, output, state -> 
        [ output: output.output ]
      },
      auto: [publish: true]
    )

    | niceView()

  emit:
  output_ch
}