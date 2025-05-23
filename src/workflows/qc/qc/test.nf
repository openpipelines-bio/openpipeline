nextflow.enable.dsl=2

include { qc } from params.rootDir + "/target/_test/nextflow/workflows/qc/qc/main.nf"
include { qc_test } from params.rootDir + "/target/nextflow/test_workflows/qc/qc_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = 
    Channel.fromList([
      [
        id: "mouse_test",
        input: resources_test.resolve("concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
        var_name_mitochondrial_genes: "mitochondrial",
        var_name_ribosomal_genes: "ribosomal",
      ],
      [
        id: "human_test",
        input: resources_test.resolve("concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
        var_name_mitochondrial_genes: "mitochondrial",
        var_name_ribosomal_genes: "ribosomal",
      ]
    ])
    | map { state -> [state.id, state] }
    | qc.run(
      toState: { id, output, state -> output + [og_input: state.input] }
    )

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
    | qc_test.run(
      fromState: {id, state ->
        [
         input: state.output,
         og_input: state.og_input
        ]
      }
    )
  
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["human_test", "mouse_test"]
    }
}
