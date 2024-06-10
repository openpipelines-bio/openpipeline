workflow process_reference {
  take:
    input_ch

  main:
    reference_ch = input_ch
    // Create reference specific output for this channel
    | map {id, state ->
      def new_state = state + ["reference_processed": state.output]
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

    | flatMap {id, state ->
        def outputDir = state.output
        def files = readCsv(state.output_files.toUriString())
        files.collect{ dat ->
          def new_id = id + "_" + dat.name
          def new_data = outputDir.resolve(dat.filename)
          [ new_id, state + ["reference_input": new_data]]
        }
        }
  
    // Remove arguments from split samples from state

    | map {id, state -> 
      def keysToRemove = ["output_files"]
      def newState = state.findAll{it.key !in keysToRemove}
      [id, newState]
    }

    | view {"After splitting samples: $it"}

    | process_samples_workflow.run(
      fromState: {id, state ->
        def newState = [
          "input": state.reference_input, 
          "id": id,]
      },
      toState: ["reference_processed": "output"]
    )


    | view {"After splitting samples: $it"}

  emit:
    reference_ch
  }

workflow process_query {
  take:
    input_ch

  main:
    query_ch = input_ch
    // Create query specific output for this channel
    | map {id, state ->
      def new_state = state + ["query_processed": state.output]
      [id, new_state]
    }

    | flatMap {id, state ->
      def outputDir = state.query_output
      def query_files = state.input
      // Workflow can take multiple input files. Split into seperate events.
      query_files.collect{ dat ->
        def filename = dat.getName()
        // make id's unique based on filename
        def new_id = id + "_" + filename.substring(0, filename.lastIndexOf('.h5mu'))
        def new_data = dat
        [ new_id, state + ["query_input": new_data]]
      }
    }
    | view {"After splitting input: $it"}

    | process_samples_workflow.run(
      fromState: {id, state ->
        def newState = [
          "input": state.query_input, 
          "id": id]
      },
      toState: ["query_processed": "output"]
    )

    | niceView()

  emit:
    query_ch
}

workflow run_wf {
  take:
    input_ch

  main:
    reference_ch = process_reference(input_ch)
    query_ch = process_query(input_ch)

    output_ch = reference_ch.mix(query_ch)
    | niceView()

  emit:
    output_ch
}