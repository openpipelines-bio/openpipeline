nextflow.enable.dsl=2

include { celltypist } from params.rootDir + "/target/nextflow/workflows/annotation/celltypist/main.nf"
include { celltypist_test } from params.rootDir + "/target/_test/nextflow/test_workflows/annotation/celltypist_test/main.nf"
params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "reference_dataset_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        reference_var_gene_names: "ensemblid",
        input_obs_batch_label: "sample_id",
        reference_obs_batch_label: "donor_assay",
        reference_obs_target: "cell_type",
      ]
    ])
    | map{ state -> [state.id, state] }
    | celltypist 
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
    | celltypist_test.run(
        fromState: [
          "input": "output"
        ],
        args: [
          "expected_modalities": ["rna", "prot"]
        ]
    )
  }

workflow test_wf_2 {
  // allow changing the resources_test dir
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "reference_model_test",
        input: resources_test.resolve("annotation_test_data/demo_2000_cells.h5mu"),
        model: resources_test.resolve("annotation_test_data/celltypist_model_Immune_All_Low.pkl"),
        reference_obs_target: "cell_type",
      ],
    ])
    | map{ state -> [state.id, state] }
    | celltypist 
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
    | celltypist_test.run(
        fromState: [
          "input": "output"
        ],
        args: [
          "expected_modalities": ["rna"]
        ]
    )
  }
