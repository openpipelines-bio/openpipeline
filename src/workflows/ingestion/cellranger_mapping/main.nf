workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | cellranger_count.run(
      fromState: [
        "input": "input",
        "output": "output_raw",
        "expect_cells": "expect_cells",
        "chemistry": "chemistry",
        "secondary_analysis": "secondary_analysis",
        "generate_bam": "generate_bam",
        "include_introns": "include_introns",
        "reference": "reference"
      ],
      toState: [
        "input": "output",
        "output_raw": "output"
      ],
      auto: [ publish: true ]
    )
    // split output dir into map
    | cellranger_count_split.run(
      fromState: {id, state -> 
        def stateMapping = [
          "input": state.input,
        ]
        outputType = state.output_type == "raw" ? "raw_h5" : "filtered_h5"
        stateMapping += [outputType: "\$id.\$key.${outputType}.h5"]
        stateMapping += ["metrics_summary": "\$id.\$key.metrics_summary.csv"]
        return stateMapping
      },
      toState: {id, output, state -> 
        def outputFile = state.output_type == "raw" ? output.raw_h5 : output.filtered_h5
        def newState = state + [ "input": outputFile ] 
        return newState
      }
    )
    // convert to h5mu
    | from_10xh5_to_h5mu.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "output_compression": "gzip",
          "output": state.output_h5mu,
          "uns_metrics": state.uns_metrics,
          "input_metrics_summary": state.metrics_summary
        ]
      },
      toState: { id, output, state ->
        [
          "output_raw": state.output_raw,
          "output_h5mu": output.output
        ]
      },
      auto: [ publish: true ],
    )

  emit:
  output_ch
}