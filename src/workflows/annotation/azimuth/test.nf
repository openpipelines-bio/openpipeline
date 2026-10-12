nextflow.enable.dsl=2

include { azimuth } from params.rootDir + "/target/nextflow/workflows/annotation/azimuth/main.nf"
include { azimuth_test } from params.rootDir + "/target/_test/nextflow/test_workflows/annotation/azimuth_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file(params.resources_test)

  // What the output of each test is checked for by azimuth_test
  def expected = [
    "h5mu_input_test": [
      "expected_modalities": ["rna", "prot", "prot_imputed"],
      "expected_obs_predictions": ["azimuth_pred_celltype_l1", "azimuth_pred_celltype_l2"],
      "expected_obs_probability": ["azimuth_probability_celltype_l1", "azimuth_probability_celltype_l2"]
    ],
    "seurat_input_test": [
      "expected_modalities": ["rna"],
      "expected_obs_predictions": ["custom_pred"],
      "expected_obs_probability": ["custom_probability"],
      "expected_obs_mapping_score": "custom_mapping_score",
      "expected_obsm": ["X_custom_umap", "X_azimuth_integrated"]
    ]
  ]

  output_ch = Channel.fromList(
    [
      // H5MU input with default output names: the untouched 'prot' modality
      // and the imputed protein expression should both end up in the output.
      [
        id: "h5mu_input_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        input_var_gene_names: "gene_symbol",
        reference: "pbmcref",
        annotation_levels: ["celltype.l1", "celltype.l2"],
        do_adt: true
      ],
      // Seurat input with custom output names and a Seurat output file.
      [
        id: "seurat_input_test",
        input_seurat: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.rds"),
        input_layer: "X",
        input_var_gene_names: "gene_symbol",
        reference: "pbmcref",
        annotation_levels: ["celltype.l1"],
        output_obs_predictions: ["custom_pred"],
        output_obs_probability: ["custom_probability"],
        output_obs_mapping_score: "custom_mapping_score",
        output_obsm_umap: "X_custom_umap",
        output_obsm_integrated: "X_azimuth_integrated",
        output_seurat: "output.rds"
      ]
    ])
    | map{ state -> [state.id, state] }
    | azimuth
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

      // --output_seurat is only returned when it was provided
      if (id == "seurat_input_test") {
        assert state.containsKey("output_seurat") : "Output should contain key 'output_seurat'."
        assert state.output_seurat.isFile() : "'output_seurat' should be a file."
        assert state.output_seurat.toString().endsWith(".rds") : "Seurat output should end with '.rds'. Found: ${state.output_seurat}"
      } else {
        assert !state.containsKey("output_seurat") : "Output should not contain key 'output_seurat'. Found: ${state.output_seurat}"
      }

      "Output: $output"
    }
    | azimuth_test.run(
      fromState: {id, state -> ["input": state.output] + expected[id]}
    )
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["h5mu_input_test", "seurat_input_test"]
    }
}
