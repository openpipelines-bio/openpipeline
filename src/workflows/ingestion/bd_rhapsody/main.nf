workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // run bd rhapsody
    | bd_rhapsody_component.run(
      fromState: { id, state ->
        // pass all arguments except:
        //  - remove output_h5mu
        //  - rename output_raw to output_dir
        def data_ = state.clone()
        data_.remove("output")
        data_.remove("output_raw")
        data_
      },
      toState: [
        "input": "output_mudata",
        "input_vdj_contigs": "vdj_unfiltered_contigs",
        "output_raw": "output_dir"
      ]
    )

    // convert to h5mu
    | from_bdrhap_to_h5mu.run(
      fromState: { id, state ->
        [
          id: id,
          input: state.input,
          input_vdj_contigs: state.input_vdj_contigs,
          output: state.output,
          output_compression: "gzip"
        ]
      },
      toState: ["output": "output"]
    )

    | setState(["output_raw", "output"])

  emit:
  output_ch
}
