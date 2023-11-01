workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // run bd rhapsody
    | bd_rhapsody_component.run(
      auto: [ publish: true ],
      fromState: { id, state ->
        // pass all arguments except:
        //  - remove output_h5mu and output_compression
        //  - rename output_raw to output
        def data_ = state.clone()
        data_.remove("output_h5mu")
        data_.remove("output_raw")
        data_.remove("output_compression")
        data_ + [ output: state.output_raw ]
      },
      toState: { id, data, state ->
        state + [ output_raw: data.output ]
      }
    )
    | view {"After bd_rhapsody: $it"}

    // convert to h5mu
    | from_bdrhap_to_h5mu.run(
      fromState: { id, state ->
        [
          id: id,
          input: state.output_raw,
          output: state.output_h5mu,
          output_compression: "gzip"
        ]
      },
      toState: { id, data, state ->
        [
          output_raw: state.output_h5mu,
          output_h5mu: data.output
        ]
      },
      auto: [publish: true]
    )

  emit:
  output_ch
}
