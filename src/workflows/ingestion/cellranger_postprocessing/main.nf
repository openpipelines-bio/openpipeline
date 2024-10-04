workflow run_wf {
  take:
  input_ch

  main:
  // perform correction if so desired

  output_ch = input_ch
    | cellbender_remove_background.run(
      runIf: {id, state -> state.perform_correction},
      fromState: { id, state ->
        [
          input: state.input,
          epochs: state.cellbender_epochs,
          output_layer: "cellbender_corrected",
          output_compression: "gzip"
        ]
      },
      toState: { id, output, state -> 
        state + ["input": output.output, "layer": "cellbender_corrected"]
      }
    )
    | filter_with_counts.run(
      runIf: {id, state -> 
        state.min_genes != null || state.min_counts != null
      },
      fromState: { id, state ->
        [
          input: state.input,
          min_genes: state.min_genes,
          min_counts: state.min_counts,
          layer: state.layer,
          output_compression: "gzip",
          do_subset: true
        ]
      },
      toState: [input: "output"]
    )
    // Make sure to use the correct ouput file names, 
    // irrespective wether or not any of the above 
    // components were run
    | publish.run(
      fromState: [ input: "input", output: "output" ],
    )
    | setState(["output": "input"])

  emit:
  output_ch
}
