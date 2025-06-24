nextflow.enable.dsl=2

include { split_h5mu } from params.rootDir + "/target/_private/nextflow/workflows/multiomics/split_h5mu/main.nf"
include { split_h5mu_test } from params.rootDir + "/target/_test/nextflow/test_workflows/multiomics/split_h5mu_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
      id: "split_test",
      input: resources_test.resolve("concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu"),
      publish_dir: "foo/",
      obs_feature: "sample_id",
      output: "samples",
      output_files: "samples.csv"
    ]
  ])
  | map { state -> [state.id, state]}
  | split_h5mu.run(
    toState: { id, output, state -> output + [orig_input: state.input] }
  )
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    "Output: $output"
  }

  | split_h5mu_test.run(
    fromState: [ 
      "input": "output_files",
      "samples_dir": "output",
      "orig_input": "orig_input",
    ]
  )
}
