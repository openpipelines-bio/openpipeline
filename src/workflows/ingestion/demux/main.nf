workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // translate the demultiplexer argument to the component name that needs to be run
    | map { id, state ->
      def funcNameMapper = [
        bclconvert: "bcl_convert",
        bcl2fastq: "bcl2fastq",
        mkfastq: "cellranger_mkfastq"
      ]
      def funcName = funcNameMapper[state.demultiplexer]
      def newState = state + [funcName: funcName]
      [id, newState]
    }

    // run the demultiplexers
    | runComponents(
      components: [cellranger_mkfastq, bcl_convert, bcl2fastq],
      filter: { id, state, config ->
        state.funcName == config.functionality.name
      },
      fromState: { id, state, config ->
        def data = [
          input: state.input,
          sample_sheet: state.sample_sheet,
          reports: null // disable reports so they end up in the output dir
        ]
        if (config.functionality.name == "bcl2fastq") {
          data.ignore_missing = state.ignore_missing
        }
        data
      },
      toState: ["output_fastq": "output"]
    )

    // run fastqc
    | fastqc.run(
      fromState: ["input": "output_fastq"],
      args: [mode: "dir"],
      toState: ["output_fastqc": "output"]
    )

    // run multiqc
    | multiqc.run(
      fromState: { id, state ->
        ["input": [state.output_fastq]]
      },
      toState: ["output_multiqc": "output"]
    )

    // subset state to the outputs
    | setState(["output_fastq", "output_fastqc", "output_multiqc"])

  emit:
  output_ch
}
