process splitStub {
  input:
    tuple val(id), path(unused)

  output:
    tuple val(id), path("stub_h5mus"), path("modalities.csv")

  script:
    """
    echo "This process is not meant to be run without -stub being defined."
    exit 1
    """

  stub:
    """
    mkdir stub_h5mus
    touch stub_h5mus/${id}_vdj.h5mu
    touch stub_h5mus/${id}_rna.h5mu
    touch stub_h5mus/${id}_prot.h5mu
    echo -e "name,filename\nrna,stub_h5mus/${id}_rna.h5mu\nprot,stub_h5mus/${id}_prot.h5mu\nvdj,stub_h5mus/${id}_vdj.h5mu" > modalities.csv
    """
}

workflow run_wf {
  // Split multimodal MuData files into several unimodal MuData files.
  take:
    input_ch

  main:
    split_ch = input_ch
      | split_modalities_component.run(
        filter: {!workflow.stubRun}, 
        fromState: ["input": "input"],
        toState: [
          "output": "output",
          "output_types": "output_types"
        ]
      )

    split_stub_ch = input_ch
      | filter{workflow.stubRun}
      // This is not a build viash component, so we cannot use
      // fromState or toState functionality
      | map {id, state -> [id, state.input]}
      | splitStub
      | map {id, output, output_types ->
        [id, ["output": output, "output_types": output_types]]
      }

    output_ch = split_ch.concat(split_stub_ch)
      | setState(["output", "output_types"])

  emit:
    output_ch
}