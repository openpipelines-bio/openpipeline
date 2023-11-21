workflowDir = params.rootDir + "/src/workflows"
include { splitStub } from workflowDir + '/multiomics/split_modalities/split_stub.nf'


workflow run_wf {
  // Split multimodal MuData files into several unimodal MuData files.
  take:
    input_ch

  main:
    split_ch = input_ch
      | split_modalities_component.run(
        runIf: {id, state -> !workflow.stubRun},
        fromState: ["input": "input"],
        toState: [
          "output": "output",
          "output_types": "output_types"
        ]
      )

    split_stub_ch = input_ch
      | filter{workflow.stubRun}
      | map {id, state -> [id, state.input, state]}
      // This is not a build viash component, so we cannot use
      // fromState or toState functionality
      | splitStub
      | map {id, output, output_types, state -> 
        [id, state + ["output": output, "output_types": output_types]]
      }

    output_ch = split_ch.concat(split_stub_ch)
      | setState(["output", "output_types"])

  emit:
    output_ch
}