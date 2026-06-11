nextflow.enable.dsl=2

include { qc_filter } from params.rootDir + "/target/nextflow/workflows/qc/qc_filter/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "qc_filter_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        rna_min_count: 100,
        rna_max_quantile: 0.99,
        max_pct_counts_mt: 50,
        var_gene_names: "gene_symbol",
        prot_min_count: 5,
        prot_max_quantile: 0.99,
        filter_genes_min_cells: 3,
        output: "qc_filter_test.output.h5mu",
        cell_count_report: "qc_filter_test.cell_counts.tsv"
      ],
    ])
    | map{ state -> [state.id, state] }
    | qc_filter
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, state]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output should be a h5mu file. Found: ${output[1].output}"
      assert output[1].cell_count_report.toString().endsWith(".tsv") : "Cell count report should be a tsv file. Found: ${output[1].cell_count_report}"
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "qc_filter_test" : "Output ID should be 'qc_filter_test'"
    }
}

workflow test_wf_skip_scrublet {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "qc_filter_skip_scrublet_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        rna_min_count: 100,
        rna_max_quantile: 0.99,
        max_pct_counts_mt: 50,
        var_gene_names: "gene_symbol",
        prot_min_count: 5,
        prot_max_quantile: 0.99,
        filter_genes_min_cells: 3,
        skip_scrublet: true,
        output: "qc_filter_skip_scrublet_test.output.h5mu",
        cell_count_report: "qc_filter_skip_scrublet_test.cell_counts.tsv"
      ],
    ])
    | map{ state -> [state.id, state] }
    | qc_filter
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, state]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output should be a h5mu file. Found: ${output[1].output}"
      assert output[1].cell_count_report.toString().endsWith(".tsv") : "Cell count report should be a tsv file. Found: ${output[1].cell_count_report}"
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "qc_filter_skip_scrublet_test" : "Output ID should be 'qc_filter_skip_scrublet_test'"
    }
}
