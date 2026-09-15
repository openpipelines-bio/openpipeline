nextflow.enable.dsl=2

include { trajectory_analysis } from params.rootDir + "/target/nextflow/workflows/beyond/trajectory_analysis/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  resources_test = file(params.resources_test)

  Channel.fromList([
    [
      id:    "beyond_atlas",
      input: resources_test.resolve("beyond_test_data/atlas.h5mu"),

      output: "beyond_atlas.h5mu",

      // Atlas metadata
      obs_participant_id: "participant_id",
      obs_subpopulation:  "subpopulation",
      obs_cluster:        "leiden",

      // Palantir: direct cell barcode as root; real-data waypoints (Palantir caps to n_cells)
      start_cell:    "cell_00000",
      num_waypoints: 500,

      // VIA: use ExN.1 as root cluster (matches leiden/subpopulation category names)
      root_user: "ExN.1",

      // Pseudotime dynamics
      obs_pseudotime: "palantir_pseudotime",
      n_splines:      4,
      dynamics_lam:   0.6,
      uns_dynamics:   "dynamics",

      // Cellular communities (small n for test speed)
      n_communities:       2,
      communities_alpha:   0.5,
      communities_method:  "hierarchical",

      // Trait associations (use continuous traits that vary across donors in the test CSV)
      traits_csv:     resources_test.resolve("beyond_test_data/traits.csv"),
      trait_columns:  ["age", "pmi"],

      // Pathway enrichment (local GMT file so no internet required in CI)
      de_results_csv:          resources_test.resolve("beyond_test_data/de_ExN.csv"),
      // gene names are the row index of the CSV (no named gene_column needed)
      gene_sets:               [resources_test.resolve("beyond_test_data/gene_sets.gmt").toString()],
      pathway_method:          "prerank",
      output_pathway_csv_dir:  "pathway_enrichment_results",
    ]
  ])
  | map { state -> [state.id, state] }
  | trajectory_analysis
  | view { output ->
    assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

    def state = output[1]
    assert state instanceof Map : "State should be a map. Found: ${state}"
    assert state.containsKey("output") : "Output should contain key 'output'."
    assert state.output.isFile() : "'output' should be a file."
    assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

    "Output: $output"
  }
}
