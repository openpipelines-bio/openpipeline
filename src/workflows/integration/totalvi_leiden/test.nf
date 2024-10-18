nextflow.enable.dsl=2

include { totalvi_leiden } from params.rootDir + "/target/nextflow/workflows/integration/totalvi_leiden/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

    output_ch = Channel.fromList([
      [
        id: "simple_execution_test",
        input: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        reference: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        prot_modality: "prot",
        prot_reference_modality: "prot",
        var_input: "filter_with_hvg",
        reference_model_path: "totalvi_reference_model",
        query_model_path: "totalvi_query_model",
        max_epochs: 1,
        max_query_epochs: 1,
      ],
      [
        id: "no_prot_leiden_resolutions_test",
        input: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        reference: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        prot_modality: "prot",
        prot_reference_modality: "prot",
        var_input: "filter_with_hvg",
        reference_model_path: "totalvi_reference_model",
        query_model_path: "totalvi_query_model",
        max_epochs: 1,
        max_query_epochs: 1,
        prot_leiden_resolution: []
      ],
      [
        id: "no_rna_leiden_resolutions_test",
        input: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        reference: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        prot_modality: "prot",
        prot_reference_modality: "prot",
        var_input: "filter_with_hvg",
        reference_model_path: "totalvi_reference_model",
        query_model_path: "totalvi_query_model",
        max_epochs: 1,
        max_query_epochs: 1,
        rna_leiden_resolution: []
      ]
    ])
    | map{ state -> [state.id, state] }
    | totalvi_leiden
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

      // check reference_model
      assert state.containsKey("reference_model_path") : "Output should contain key 'reference_model_path'."
      assert state.reference_model_path.isDirectory() : "'reference_model_path' should be a directory."
      assert state.reference_model_path.toString().endsWith("_reference_model") : "Model output directory should end with '_reference_model'. Found: ${state.reference_model_path}"

      // check query_model
      assert state.containsKey("query_model_path") : "Output should contain key 'query_model_path'."
      assert state.query_model_path.isDirectory() : "'query_model_path' should be a directory."
      assert state.query_model_path.toString().endsWith("_query_model") : "Model output directory should end with '_query_model'. Found: ${state.query_model_path}" 
      
      
      "Output: $output"
    }
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 3 : "output channel should contain 3 events"
      assert output_list.collect{it[0]} == ["no_prot_leiden_resolutions_test", "no_rna_leiden_resolutions_test", "simple_execution_test"]
    }
}
