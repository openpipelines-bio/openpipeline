workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map{id, state ->
      def output_state
      if (!workflow.stubRun) {
        assert (state.perform_correction || state.min_genes != null || state.min_counts != null):
          "Either perform_correct, min_genes or min_counts should be specified!"
        output_state = state
      } else {
        output_state = state + ["perform_correction": true, "min_genes": 1, "min_counts": 1]
      }
      [id, output_state]
    }
    // Make sure there is not conflict between the output from this workflow
    // And the output from any of the components
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // perform correction if so desired
    | cellbender_remove_background.run(
      runIf: {id, state -> state.perform_correction},
      fromState: { id, state ->
        [
          input: state.input,
          epochs: state.cellbender_epochs,
          output_layer: "cellbender_corrected",
          output_compression: "gzip",
          output: state.workflow_output,
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
          do_subset: true,
          output: state.workflow_output,
        ]
      },
      toState: [input: "output"]
    )
    // Make sure to use the correct ouput file names, 
    // irrespective of which component(s) (or combinations of then)
    // were run. The above components
    // should put their output into 'input'
    | map {id, state -> 
      [id, ["output": state.input]]
    }
    | view {"HERE: $it"}

  emit:
  output_ch
}
