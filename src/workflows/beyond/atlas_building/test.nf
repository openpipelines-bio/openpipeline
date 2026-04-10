nextflow.enable.dsl=2

include { atlas_building } from params.rootDir + "/target/nextflow/workflows/beyond/atlas_building/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  resources_test = file(params.resources_test)

  shared = [
    output:                                  "atlas.h5mu",
    reference:                               resources_test.resolve("beyond_test_data/reference.h5mu"),
    reference_obs_target:                    "cell_type",
    obs_covariates:                          ["participant_id"],
    highly_variable_features_obs_batch_key:  "participant_id",
    leiden_resolution:                       [0.5],
    subpop_leiden_resolution:                [0.3],
    skip_scrublet_doublet_detection:         true,
    min_counts:                              1,
    min_genes_per_cell:                      1,
  ]

  Channel.fromList([
    [ id: "donor_01", input: resources_test.resolve("beyond_test_data/donor_01.h5mu") ] + shared,
    [ id: "donor_02", input: resources_test.resolve("beyond_test_data/donor_02.h5mu") ] + shared,
    [ id: "donor_03", input: resources_test.resolve("beyond_test_data/donor_03.h5mu") ] + shared,
    [ id: "donor_04", input: resources_test.resolve("beyond_test_data/donor_04.h5mu") ] + shared,
  ])
  | map { state -> [state.id, state] }
  | atlas_building
  | view { output ->
    assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

    def state = output[1]
    assert state instanceof Map : "State should be a map. Found: ${state}"
    assert state.containsKey("output") : "Output should contain key 'output'."
    assert state.output.isFile() : "'output' should be a file."
    assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

    "Output: $output"
  }
}
