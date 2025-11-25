nextflow.enable.dsl=2

include { scanvi_model } from params.rootDir + "/target/nextflow/workflows/embedding_model/scanvi_model/main.nf"
include { scanvi_model_test } from params.rootDir + "/target/_test/nextflow/test_workflows/embedding_model/scanvi_model_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "simple_execution_test",
        input: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        obs_batch_label: "donor_assay",
        obs_target: "cell_type",
        var_gene_names: "ensemblid",
        max_epochs: 10
      ]
    ])
    | map{ state -> [state.id, state] }
    | scanvi_model 
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test") : "Output ID should be same as input ID"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.containsKey("output_scanvi_model") : "Output should contain key 'output_scanvi_model'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

      // check output_model
      assert state.containsKey("output_scanvi_model") : "Output should contain key 'output_scanvi_model'."
      assert state.output_scanvi_model.isDirectory() : "'output_scanvi_model' should be a directory."
      assert state.output_scanvi_model.toString().endsWith("_model") : "Model output directory should end with '_model'. Found: ${state.output_scanvi_model}"
      def modelFile = state.output_scanvi_model.resolve("model.pt")
      assert modelFile.isFile(), "Model output directory should contain a model.pt file."

    "Output: $output"
    }
    | scanvi_model_test.run(
        fromState: [
          "input": "output"
        ]
    )

  output_ch
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert output_list.collect{it[0]} == ["simple_execution_test"]
    }
  }

  workflow test_wf_2 {
  // allow changing the resources_test dir
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "simple_execution_test",
        obs_covariate: ["assay"],
        input: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        obs_batch_label: "donor_assay",
        obs_target: "cell_type",
        var_gene_names: "ensemblid",
        max_epochs: 10
      ]
    ])
    | map{ state -> [state.id, state] }
    | scanvi_model 
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test") : "Output ID should be same as input ID"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.containsKey("output_scanvi_model") : "Output should contain key 'output_scanvi_model'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

      // check output_model
      assert state.containsKey("output_scanvi_model") : "Output should contain key 'output_scanvi_model'."
      assert state.output_scanvi_model.isDirectory() : "'output_scanvi_model' should be a directory."
      assert state.output_scanvi_model.toString().endsWith("_model") : "Model output directory should end with '_model'. Found: ${state.output_scanvi_model}"
      def modelFile = state.output_scanvi_model.resolve("model.pt")
      assert modelFile.isFile(), "Model output directory should contain a model.pt file."

    "Output: $output"
    }
    | scanvi_model_test.run(
        fromState: [
          "input": "output"
        ]
    )

  output_ch
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert output_list.collect{it[0]} == ["simple_execution_test"]
    }
  }