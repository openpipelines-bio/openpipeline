workflow run_wf {
  take:
  input_ch

  main:
  // perform correction if so desired
  mid1_corrected = input_ch
    | filter{ it[1].perform_correction }
    | cellbender_remove_background.run(
      fromState: { id, state ->
        [
          input: state.input,
          epochs: state.cellbender_epochs,
          output_layer: "cellbender_corrected",
          output_compression: "gzip"
        ]
      },
      toState: { id, output, state -> 
        state + [input: output.output, layer: "cellbender_corrected"]
      }
    )
  mid1_uncorrected = input_ch
    | filter{ ! it[1].perform_correction }
  mid1 = mid1_corrected.mix(mid1_uncorrected)

  // perform filtering if so desired
  mid2_filtered = mid1
    | filter{ it[1].min_genes != null || it[1].min_counts != null }
    | filter_with_counts.run(
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
  mid2_unfiltered = mid1
    | filter{ it[1].min_genes == null && it[1].min_counts == null }
  mid2 = mid2_filtered.mix(mid2_unfiltered)
    
  // return output map
  output_ch = mid2
    | publish.run(
      fromState: [ input: "input", output: "output" ],
      auto: [ publish: true ]
    )

  emit:
  output_ch
}
