nextflow.enable.dsl=2

include { scgpt_integration_knn } from params.rootDir + "/target/nextflow/workflows/annotation/scgpt_integration_knn/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList(
    [
      [
        id: "simple_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms_w_sample_id.h5mu"),
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        model: resources_test.resolve("scgpt/source/best_model.pt"),
        model_config: resources_test.resolve("scgpt/source/args.json"),
        model_vocab: resources_test.resolve("scgpt/source/vocab.json"),
        obs_reference_targets: "cell_type",
        obs_batch_label: "donor_assay",
        var_gene_names_query: "gene_symbol",
        n_hvg: 400,
        seed: 1,
        leiden_resolution: [1.0, 0.25]
      ],
      [
        id: "no_leiden_resolutions_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms_w_sample_id.h5mu"),
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        model: resources_test.resolve("scgpt/source/best_model.pt"),
        model_config: resources_test.resolve("scgpt/source/args.json"),
        model_vocab: resources_test.resolve("scgpt/source/vocab.json"),
        var_gene_names_query: "gene_symbol",
        obs_reference_targets: "cell_type",
        obs_batch_label: "donor_assay",
        n_hvg: 400,
        seed: 1,
        leiden_resolution: []
      ]
    ])
    | map{ state -> [state.id, state] }
    | scgpt_integration_knn 
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
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

    "Output: $output"
    }
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["no_leiden_resolutions_test", "simple_execution_test"]
    }

    }