


workflow run_wf {
    take: 
    input_ch
        
    main:
    output_ch = input_ch
      | normalize_total.run(
        fromState: [
          "input": "input",
          "modality": "modality",
          "input_layer": "layer",
          "target_sum": "target_sum"
        ],
        args: [
          "output_layer": "normalized",
        ],
        toState: [
          "input": "output",
        ]
      )
      | log1p.run( 
        fromState: [
          "input": "input",
          "modality": "modality",
          "output_layer": "output_layer"
        ],
        args: [
          "input_layer": "normalized",
        ],
        toState: [
          "input": "output"
        ]
      )
      | delete_layer.run(
        fromState: [
          "input": "input",
          "modality": "modality"
        ],
        args: [
          "layer": "normalized"
        ],
        toState: [
          "output": "output"
        ]
      )
      | setState(["output"])

    emit:
    output_ch


}