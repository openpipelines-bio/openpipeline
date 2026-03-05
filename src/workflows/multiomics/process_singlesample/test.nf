nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { process_singlesample } from targetDir + "/workflows/multiomics/process_singlesample/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  input_ch = Channel.fromList([
    [
      id: "mouse",
      input: resources_test.resolve("concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      publish_dir: "foo/",
      output: "test.h5mu",
    ],
    [
      id: "human",
      input: resources_test.resolve("concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      publish_dir: "foo/",
      output: "test.h5mu",

    ]
  ])
  | map{ state -> [state.id, state] }
  | process_singlesample.run(
    toState: { id, output, state -> output + [orig_input: state.input] }
  )

  assert_ch = input_ch
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith("test.h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 2 : "output channel should contain two events, got ${output_list.size()}"
    assert output_list.collect({it[0]}).sort() == ["human", "mouse"] : "Output id's should contain `mouse` and `human`.'"
  }
      
  // test_ch = input_ch
  //   | workflow_test.run(
  //     fromState: [
  //       "input": "output",
  //       "orig_input": "orig_input"
  //     ],
  //   )
  
}