nextflow.enable.dsl=2

include { conversion } from params.rootDir + "/target/_test/nextflow/workflows/ingestion/conversion/main.nf"
include { conversion_test } from params.rootDir + "/target/nextflow/test_workflows/ingestion/conversion_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "10xh5_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5"),
        input_type: "10xh5",
        modality: null
      ],
      [
        id: "10xmtx_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix"),
        input_type: "10xmtx",
        modality: null,
        output: "\$id.h5mu"
      ],
      [
        id: "10xmtx",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix"),
        input_type: "10xmtx",
        modality: "rna",
        output: "\$key.h5mu"
      ],
      [
        id: "h5ad",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5ad"),
        input_type: "h5ad",
        modality: "rna",
        output: "\$key.h5mu"
      ]
    ])
    | map{ state -> [state.id, state] }
    | conversion 
    | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
    }

    | conversion_test.run(
      fromState: ["input": "output"]
    )

    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 4 : "output channel should contain four events"
    }
}
