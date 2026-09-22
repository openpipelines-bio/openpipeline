nextflow.enable.dsl=2

include { subpopulation_clustering } from params.rootDir + "/target/nextflow/workflows/beyond/subpopulation_clustering/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  resources_test = file(params.resources_test)

  Channel.fromList([
    [
      id:             "beyond_atlas",
      input:          resources_test.resolve("beyond_test_data/atlas.h5mu"),
      output:         "atlas_with_subpopulations.h5mu",
      obs_cell_type:  "celltypist_pred",
      obs_label:      "leiden_subpopulation",
      leiden_resolution: [0.3],
    ]
  ])
  | map { state -> [state.id, state] }
  | subpopulation_clustering
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
