nextflow.enable.dsl=2

include { pseudobulk_deseq2 } from params.rootDir + "/target/nextflow/workflows/differential_expression/pseudobulk_deseq2/main.nf"
params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "simple_execution_test",
        input: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        obs_label: "cell_type",
        obs_groups: "treatment",
        design_formula: "~ cell_type + treatment",
        contrast_column: "treatment",
        contrast_values: ["ctrl", "stim"],
        output: "simple_execution_test_output.csv"
      ]
    ])
    | map{ state -> [state.id, state] }
    | pseudobulk_deseq2
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test") : "Output ID should be same as input ID"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".csv") : "Output file should end with '.csv'. Found: ${state.output}"
    
    "Output: $output"
    }
}