workflow run_wf {
  
  take:
    input_ch

  main:
    output_ch = input_ch
    // Set aside the output for this workflow to avoid conflicts
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }

    | binning.run(
        fromState: {id, state -> [
            "input": state.input,
            "output": "workflow_output"
          ]
        },
        auto: [ publish: true ]
    )
    
  emit:
    output_ch
}
