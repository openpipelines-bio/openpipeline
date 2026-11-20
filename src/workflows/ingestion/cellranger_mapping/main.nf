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
        "disable_cell_annotation": "disable_cell_annotation",
        "secondary_analysis": "secondary_analysis",
        "generate_bam": "generate_bam",
        "include_introns": "include_introns",
        "reference": "reference"
      ],
      toState: ["output_raw": "output"]
    )
    // convert to h5mu
    | from_10xh5_to_h5mu.run(
      fromState: {id, state ->
        def h5File = state.output_type == "raw" ? "raw_feature_bc_matrix.h5" : "filtered_feature_bc_matrix.h5"
        [
          "input": state.output_raw.resolve(h5File),
          "output_compression": "gzip",
          "output": state.output_h5mu,
          "uns_metrics": state.uns_metrics,
          "input_metrics_summary": state.output_raw.resolve("metrics_summary.csv")
        ]
      },
      toState: ["output_h5mu": "output"]
    )
    | setState(["output_raw", "output_h5mu"])

  emit:
  output_ch
}