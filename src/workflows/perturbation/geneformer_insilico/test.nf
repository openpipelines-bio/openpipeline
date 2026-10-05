nextflow.enable.dsl=2

include { geneformer_leiden } from params.rootDir + "/target/nextflow/workflows/integration/geneformer_leiden/main.nf"
include { geneformer_insilico } from params.rootDir + "/target/nextflow/workflows/perturbation/geneformer_insilico/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)
  model = resources_test.resolve("geneformer/Geneformer-V1-10M")

  // the perturbation starts from the embedded and clustered output of
  // geneformer_leiden, so run that first
  output_ch = Channel.fromList([
      [
        id: "simple_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        model: model,
        model_version: "V1",
        leiden_resolution: [1.0],
      ],
    ])
    | map{ state -> [state.id, state] }
    | geneformer_leiden
    | map { id, state ->
      [id, [
        input: state.output,
        model: model,
        model_version: "V1",
        var_gene_names: "gene_symbol",
        obs_cluster: "geneformer_integration_leiden_1.0",
        disease_clusters: ["0"],
        healthy_clusters: ["1"],
        skip_qc: true,
        min_genes: 300,
        n_disease_cells: 4,
        batch_size: 2,
        top_n_genes: 5,
        min_coverage: 1,
        n_random: 10,
        top_n: 3,
      ]]
    }
    | geneformer_insilico
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id == "simple_execution_test" : "the batches should be gathered back into one event, got $id"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu")

      // 4 cells in batches of 2
      assert state.output_shifts instanceof List : "'output_shifts' should be a list"
      assert state.output_shifts.size() == 2 : "expected 2 batches, got ${state.output_shifts.size()}"
      def n_knockouts = state.output_shifts.sum { it.readLines().size() - 1 }
      assert n_knockouts == 4 * 5 : "expected 5 knockouts for each of 4 cells, got $n_knockouts"

      def ranked = state.output_ranked.readLines()
      assert ranked[0] == "rank,gene_id,gene_name,median_shift,pvalue,n_cells"
      assert ranked.size() > 1 : "no gene was ranked"
      def top = state.output_top.readLines()
      assert top.size() == Math.min(4, ranked.size()) : "--top_n 3 should give 3 genes"

      "Output: $output"
    }
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
    }
}
