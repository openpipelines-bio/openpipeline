process splitStub {
  input:
    tuple val(id), val(unused)

  output:
    tuple val(id), path("stub_h5mus"), path("samples.csv")

  script:
    """
    echo "This process is not meant to be run without -stub being defined."
    exit 1
    """

  stub:
    """
    mkdir stub_h5mus
    touch stub_h5mus/${id}_sample_1.h5mu
    touch stub_h5mus/${id}_sample_2.h5mu
    touch stub_h5mus/${id}_sample_3.h5mu
    echo -e "name,filename\nsample_1,stub_h5mus/${id}_sample_1.h5mu\nsample_2,stub_h5mus/${id}_sample_2.h5mu\nsample_3,stub_h5mus/${id}_sample_3.h5mu" > samples.csv
    """
}

workflow run_wf {
  take:
    input_ch

  main:
    split_ch = input_ch
      | split_h5mu_component.run(
        filter: {!workflow.stubRun}, 
        fromState: [
            "input": "input",
            "output_compression": "output_compression",
            "output_files": "output_files",
            "modality": "modality",
            "obs_feature": "obs_feature",
            "drop_obs_nan": "drop_obs_nan",
            "ensure_unique_filenames": "ensure_unique_filenames"
        ],
        toState: [
          "output": "output",
          "output_files": "output_files"
        ]
      )

    split_stub_ch = input_ch
      | filter{workflow.stubRun}
      // This is not a build viash component, so we cannot use
      // fromState or toState functionality
      | map {id, state -> [id, state.input]}
      | splitStub
      | map {id, output, output_files ->
        [id, ["output": output, "output_files": output_files]]
      }

    output_ch = split_ch.concat(split_stub_ch)
      | setState(["output", "output_files"])

  emit:
    output_ch
}