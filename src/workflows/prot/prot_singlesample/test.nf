nextflow.enable.dsl=2

include { prot_singlesample } from params.rootDir + "/target/nextflow/workflows/prot/prot_singlesample/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  output_ch = Channel.fromList([
      [
        id: "foo",
        input: file(params.resources_test).resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        min_counts: 3,
        max_counts: 100000,
        min_genes_per_cell: 2,
        max_genes_per_cell: 10000,
        min_cells_per_gene: 10,
        min_fraction_mito: 0.2,
        max_fraction_mito: 0.8,
        output: "foo.final.h5mu",
      ]
    ])
    | map{ state -> [state.id, state] }
    | prot_singlesample
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["foo.final.h5mu"] as Set)
    }
}
