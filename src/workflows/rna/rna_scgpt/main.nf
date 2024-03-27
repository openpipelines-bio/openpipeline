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

    | cross_check_genes.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": state.modality,
          "vocab_file": state.vocab_file,
          "gene_name_layer": state.gene_name_layer,
          "output": state.output,
          "pad_token": state.pad_token
        ]
      },
      toState: ["input": "output"]
    )
    | binning.run(
        fromState: {id, state -> [
            "input": state.input,
            "modality": state.modality,
            "input_layer": state.input_layer,
            "n_input_bins": state.n_input_bins,
            "output_compression": state.output_compression,
            "binned_layer": state.binned_layer,
            "output": "workflow_output"
          ]
        },
        auto: [ publish: true ]
    )
    
  emit:
    output_ch
}
