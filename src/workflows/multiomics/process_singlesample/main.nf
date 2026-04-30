workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      // Make sure there is not conflict between the output from this workflow
      // And the output from any of the components
      | map {id, state ->
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
      }
      | process_singlesample_base.run(
        fromState: {id, state -> state},
        toState: [
          "input": "output",
          "modality": "output_modality"
        ]
      )
      | groupTuple(by: 0, sort: "hash")
      | map { id, states ->
          def new_input = states.collect{it.input}
          def modalities = states.collect{it.modality}.unique()
          def other_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
          }.minus(["output", "input", "modality"])
          def new_state = other_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "Arguments should be the same across modalities. Please report this \
                                                 as a bug. Argument name: $argument_name, \
                                                 argument value: $argument_values"
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
          }
          [id, new_state + ["input": new_input, "modalities": modalities]]
      }
      | view {"Input merge channel: $it"}
      | merge.run(
        fromState: [
          "input": "input",
          "output": "workflow_output"
        ],
        toState: ["output": "output"]
      )
      | view {"After smerging: $it"}
      | intersect_obs.run(
        runIf: {id, state -> state.intersect_obs},
        fromState:[
          "input": "output",
          "modalities": "modalities",
          "output": "workflow_output"
        ],
        toState: {id, output, state ->
          ["output": output.output]
        }
      )
      | view {"After singlesample processing: $it"}
      | setState(["output"])
    // output_ch = modalities_ch

  emit:
    output_ch
}