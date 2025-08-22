nextflow.enable.dsl=2

include { scanvi_scarches } from params.rootDir + "/target/nextflow/workflows/annotation/scanvi_scarches/main.nf"
include { scanvi_scarches_test } from params.rootDir + "/target/_test/nextflow/test_workflows/annotation/scanvi_scarches_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "simple_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        input_obs_batch_label: "sample_id",
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        reference_obs_batch_label: "donor_assay",
        reference_obs_target: "cell_type",
        reference_var_gene_names: "ensemblid",
        max_epochs: 10,
        obsm_embedding: "X_pca"
      ]
    ])
    | map{ state -> [state.id, state] }
    | scanvi_scarches 
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test") : "Output ID should be same as input ID"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.containsKey("output_model") : "Output should contain key 'output_model'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

      // check output_model
      assert state.containsKey("output_model") : "Output should contain key 'output_model'."
      assert state.output_model.isDirectory() : "'output_model' should be a directory."
      assert state.output_model.toString().endsWith("_model") : "Model output directory should end with '_model'. Found: ${state.output_model}"
      def modelFile = state.output_model.resolve("model.pt")
      assert modelFile.isFile(), "Model output directory should contain a model.pt file."

    "Output: $output"
    }
    | scanvi_scarches_test.run(
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