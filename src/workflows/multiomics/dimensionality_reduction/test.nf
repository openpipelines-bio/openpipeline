nextflow.enable.dsl=2

include { dimensionality_reduction } from params.rootDir + "/target/nextflow/workflows/multiomics/dimensionality_reduction/main.nf"
include { dimensionality_reduction_test } from params.rootDir + "/target/nextflow/test_workflows/multiomics/dimensionality_reduction_test/main.nf" 

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
      [
        id: "simple_execution_test",
        input: resources_test.resolve("concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu"),
        layer: "",
        output: "foo.final.h5mu"
      ],
      [
        id: "pca_obsm_output_test",
        input: resources_test.resolve("concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu"),
        layer: "",
        output: "foo.final.h5mu"
      ],
    ])
    | map{ state -> [state.id, state] }
    | dimensionality_reduction

    assert_ch = output_ch
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test")

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

      "Output: $output"
    }

        
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["pca_obsm_output_test", "simple_execution_test"]
      output_list
    }

    test_ch = output_ch
    | map { id , output -> [id, ["input": output.output]]}
    | dimensionality_reduction_test


}