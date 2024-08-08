workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // run bd rhapsody
    | bd_rhapsody2_component.run(
      // auto: [ publish: true ],
      fromState: { id, state ->
        // pass all arguments except:
        //  - remove output_h5mu
        //  - rename output_raw to output_dir
        def data_ = state.clone()
        data_.remove("output")
        data_.remove("output_raw")
        data_ + [ output_dir: state.output_raw ]
        data_ + [ output_mudata: "bd_rhap_output.h5mu" ]
      },
      toState: [
        "input": "output_mudata",
        "output_raw": "output_dir"
      ]
    )
    | view {"After bd_rhapsody: $it"}

    // convert to h5mu
    | from_bdrhap2_to_h5mu.run(
      fromState: { id, state ->
        [
          id: id,
          input: state.input,
          output: state.output,
          output_compression: "gzip"
        ]
      },
      toState: { id, data, state ->
        [
          output_raw: state.output_raw,
          output: data.output
        ]
      },
      auto: [publish: true]
    )
    | view {"Before publish: $it"}

  emit:
  output_ch
}
