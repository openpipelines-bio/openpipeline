workflow run_wf {
  take:
  input_ch

  main:

  def targetMapping = [
    "build_cellranger_reference": "cellranger",
    "build_bdrhap_reference": "bd_rhapsody",
    "star_build_reference": "star" 
  ]

  output_ch = input_ch
    // split params for downstream components
    | make_reference_component.run(
      fromState: [
        "input": "input",
        "genome_fasta": "genome_fasta", 
        "transcriptome_gtf": "transcriptome_gtf",
        "ercc": "ercc",
        "output_fasta": "output_fasta",
        "output_gtf": "output_gtf",
        "subset_regex": "subset_regex"
      ],
      toState: [
        "output_fasta": "output_fasta",
        "output_gtf": "output_gtf"
      ],
      auto: [publish: true]
    )
    | runEach(
      components: [ 
        build_cellranger_reference,
        build_bdrhap_reference,
        star_build_reference
      ],
      filter: { id, state, component ->
        state.target.contains(targetMapping.get(component.config.functionality.name))
      },
      fromState: { id, state, component ->
        def target = targetMapping.get(component.config.functionality.name)
        def passed_state = [
          input: state.input,
          output: state.get("output_" + target),
          target: state.target,
          genome_fasta: state.output_fasta,
          transcriptome_gtf: state.output_gtf
        ]
        passed_state
      },
      toState: {id, output, state, component ->
        def target = targetMapping.get(component.config.functionality.name)
        def newState = state + ["output_$target": output.output]
        return newState
      },
      auto: [ publish: true ],
    )
    | setState(targetMapping.values().collect{"output_$it"} + ["output_fasta", "output_gtf"])
    | groupTuple(by: 0, sort: "hash")
    | map {id, state_list -> [id, state_list.inject{acc, val -> acc + val}]}
  
  emit:
  output_ch
}