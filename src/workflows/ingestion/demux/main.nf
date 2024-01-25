workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // run the demultiplexers
    | runEach(
      components: [cellranger_mkfastq, bcl_convert, bcl2fastq],
      filter: { id, state, component ->
        def funcNameMapper = [
          "bclconvert": "bcl_convert",
          "bcl2fastq": "bcl2fastq",
          "mkfastq": "cellranger_mkfastq"
        ]
        funcNameMapper[state.demultiplexer] == component.config.functionality.name
      },
      fromState: { id, state, component ->
        def data = [
          input: state.input,
          sample_sheet: state.sample_sheet,
          reports: null // disable reports so they end up in the output dir
        ]
        if (component.config.functionality.name== "bcl2fastq") {
          data.ignore_missing = state.ignore_missing
        }
        data
      },
      toState: [
        "input": "output",
        "output_fastq": "output"
      ]
    )

    // run fastqc
    | fastqc.run(
      fromState: [
        "input": "input",
        "output": "output_fastqc"
      ],
      args: [mode: "dir"],
      toState: [
        "output_fastqc": "output",
        "input": "output"
      ]
    )

    // run multiqc
    | multiqc.run(
      fromState: { id, state ->
        [
          "input": [state.input],
          "output": state.output_multiqc
        ]
      },
      toState: ["output_multiqc": "output"]
    )
    // subset state to the outputs
    | setState(["output_fastq", "output_fastqc", "output_multiqc"])


  emit:
  output_ch
}