workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | cellranger_multi_component.run(
      fromState: [
        "input": "input",
        "gex_input": "gex_input",
        "abc_input": "abc_input",
        "cgc_input": "cgc_input",
        "mux_input": "mux_input",
        "vdj_input": "vdj_input",
        "vdj_t_input": "vdj_t_input",
        "vdj_t_gd_input": "vdj_t_gd_input",
        "vdj_b_input": "vdj_b_input",
        "agc_input": "agc_input",
        "output": "output_raw",
        "cell_multiplex_sample_id": "cell_multiplex_sample_id",
        "cell_multiplex_oligo_ids": "cell_multiplex_oligo_ids",
        "cell_multiplex_description": "cell_multiplex_description",
        "gex_reference": "gex_reference",
        "feature_reference": "feature_reference",
        "vdj_reference": "vdj_reference",
        "gex_expect_cells": "gex_expect_cells",
        "gex_chemistry": "gex_chemistry",
        "gex_secondary_analysis": "gex_secondary_analysis",
        "gex_generate_bam": "gex_generate_bam",
        "gex_include_introns": "gex_include_introns",
        "library_id": "library_id",
        "library_type": "library_type",
        "library_subsample": "library_subsample",
        "library_lanes": "library_lanes",
        "vdj_inner_enrichment_primers": "vdj_inner_enrichment_primers"
      ],
      toState: [
        "output_raw": "output", 
        "input": "output"
      ],
      auto: [ publish: true ]
    )
    | from_cellranger_multi_to_h5mu.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "output": state.output_h5mu,
          "uns_metrics": state.uns_metrics,
          "output_compression": "gzip"
        ]
      },
      toState: {id, output, state -> 
        [
          "output_h5mu": output.output,
          "output_raw": state.output_raw
        ]
      },
      auto: [ publish: true ],
    )

  emit:
  output_ch
}