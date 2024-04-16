nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { process_batches } from targetDir + "/workflows/multiomics/process_batches/main.nf"
include { workflow_test } from targetDir + "/test_workflows/multiomics/process_batches/workflow_test/main.nf"
include { workflow_test2 } from targetDir + "/test_workflows/multiomics/process_batches/workflow_test2/main.nf"

workflow test_wf {

  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")
  
  input_ch = Channel.fromList([
      [
          id: "test",
          input: resources_test.resolve("concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu"),
          publish_dir: "foo/",
      ],
      [
          id: "test2",
          input: resources_test.resolve("concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu"),
          publish_dir: "foo/"
      ]
    ])
    | map{ state -> [state.id, state] }
    | view { "Input: $it" }
    | process_batches.run(
      toState: { id, output, state -> output + [orig_input: state.input] }
    )

  assert_ch = input_ch
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      print "output_list: $output_list"
      assert output_list.size() == 2 : "output channel should contain two events"
      assert output_list.collect({it[0]}).sort() == ["test", "test2"] : "First output ID should be 'test'"
    }

  test_ch = input_ch
    | workflow_test.run(
      fromState: {id, state ->
        [ input: state.output, orig_input: state.orig_input]}
      )
  
}

workflow test_wf2 {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  input_ch = Channel.fromList([
      [
          input: resources_test.resolve("10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_mms.h5mu"),
          pca_overwrite: true,
          id: "test",
          publish_dir: "foo/",
          output: "test.h5mu"
      ]
    ])
    | map{ state -> [state.id, state] }
    | view { "Input: $it" }
    | process_batches.run(
      toState: { id, output, state -> output + [orig_input: state.input] }
    )

    assert_ch = input_ch
    | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file], was $output"
    assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
    "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
    print "output_list: $output_list"
    assert output_list.size() == 1 : "output channel should contain two events"
    assert output_list.collect({it[0]}).sort() == ["test"] : "First output ID should be 'test'"
    }

    test_ch = input_ch
    | workflow_test2.run(
    fromState: {id, state ->
        [ input: state.output, orig_input: state.orig_input]}
    )

  
}