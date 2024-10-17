nextflow.enable.dsl=2

include { scvi_leiden } from params.rootDir + "/target/nextflow/workflows/integration/scvi_leiden/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  output_ch = Channel.fromList([
      [
        id: "simple_execution_test",
        input: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        layer: "log_normalized",
        obs_batch: "sample_id",
        max_epochs: 1,
        output_model: "simple_execution_test_model/"
      ],
      [
        id: "no_leiden_resolutions_test",
        input: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        layer: "log_normalized",
        obs_batch: "sample_id",
        output_model: "no_leiden_resolutions_test_model/",
        max_epochs: 1,
        leiden_resolution: []
      ] 
    ])
    | map{ state -> [state.id, state] }
    | scvi_leiden
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

      // check model_output
      assert state.containsKey("output_model") : "Output should contain key 'output_model'."
      assert state.output_model.isDirectory() : "'output_model' should be a directory."
      assert state.output_model.toString().endsWith("_model") : "Model output directory should end with '_model'. Found: ${state.output_model}"

      "Output: $output"
    }
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["no_leiden_resolutions_test", "simple_execution_test"]
    }
}

