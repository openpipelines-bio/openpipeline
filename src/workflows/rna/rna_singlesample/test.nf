nextflow.enable.dsl=2

include { rna_singlesample } from params.rootDir + "/target/nextflow/workflows/rna/rna_singlesample/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
      [
        id: "mitochondrial_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        min_counts: 3,
        max_counts: 10000000,
        min_genes_per_cell: 2,
        max_genes_per_cell: 10000000,
        min_cells_per_gene: 2,
        min_fraction_mito: 0.05,
        max_fraction_mito: 0.2,
        var_gene_names: "gene_symbol",
        var_name_mitochondrial_genes: "mitochondrial",
        obs_name_mitochondrial_fraction: "fraction_mitochondrial",
        output: "mitochondrial_test.final.h5mu"
      ],
      [
        id: "simple_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        min_counts: 3,
        max_counts: 10000000,
        min_genes_per_cell: 2,
        max_genes_per_cell: 10000000,
        min_cells_per_gene: 2,
        output: "simple_execution_test.final.h5mu"
      ]
    ])
    | map{ state -> [state.id, state] }
    | rna_singlesample
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toSortedList{a, b -> a[0] <=> b[0]}
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain two events"
      println "output_list: $output_list"
      assert output_list.collect{it[0]} == ["mitochondrial_test", "simple_execution_test"] : "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["mitochondrial_test.final.h5mu", "simple_execution_test.final.h5mu"] as Set)

    }
}

workflow test_wf2 {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
      [
        id: "test_different_fraction_column",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        min_counts: 3,
        max_counts: 10000000,
        min_genes_per_cell: 2,
        max_genes_per_cell: 10000000,
        min_cells_per_gene: 2,
        min_fraction_mito: 0.05,
        max_fraction_mito: 0.2,
        var_gene_names: "gene_symbol",
        var_name_mitochondrial_genes: "mitochondrial",
        obs_name_mitochondrial_fraction: "foobar",
        output: "mitochondrial_test.final.h5mu"
      ],
    ])
    | map{ state -> [state.id, state] }
    | rna_singlesample
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toSortedList{a, b -> a[0] <=> b[0]}
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      println "output_list: $output_list"
      assert output_list.collect{it[0]} == ["test_different_fraction_column"] : "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["mitochondrial_test.final.h5mu"] as Set)

    }
}