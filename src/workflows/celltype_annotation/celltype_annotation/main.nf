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
    // // download the reference h5ad file
    // | download_file.run(
    //   fromState: { id, state ->
    //     [
    //       "input": state.reference_url,
    //       "output": "reference.h5ad",
    //       "verbose": "true",
    //     ]
    //   },
    //   toState: [
    //     "input": "output",
    //   ]
    // )
    // // convert the reference h5ad file to h5mu
    // | from_h5ad_to_h5mu.run(
    //     fromState: { id, state ->
    //     [
    //       "input": state.input,
    //       "modality": "rna",
    //     ]
    //   },
    //   toState: [
    //     "input": "output",
    //   ]
    // )
    | split_samples.run(
        fromState: { id, state ->
        [
          "input": state.reference,
          "modality": "rna",
          "obs_feature": state.obs_reference_batch
        ]
      },
      toState:  ["output": "output", "output_files": "output_files"]
    )

    | map {id, state ->
        def outputDir = state.output
        def files = readCsv(state.output_files.toUriString())
        def new_data = files.collect{ dat -> outputDir.resolve(dat.filename)}
        def new_id = id 
        [ new_id, state + ["input": new_data]]
        }
    
    // Remove arguments from split samples from state

    | niceView()

    | map {id, state -> 
      def keysToRemove = ["output_files"]
      def newState = state.findAll{it.key !in keysToRemove}
      [id, newState]
    }

    | niceView()

  emit:
  output_ch
}