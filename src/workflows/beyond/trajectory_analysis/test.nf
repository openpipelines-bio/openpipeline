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
      obs_group: "participant_id",
      obs_label:  "subpopulation",

      // Palantir: direct cell barcode as root; real-data waypoints (Palantir caps to n_cells)
      start_cell:    "cell_00000",
      num_waypoints: 500,
      // Every cell of a donor carries that donor's proportion vector, so the PHATE
      // embedding holds one distinct point per donor, repeated once per cell. Palantir's
      // diffusion maps go singular when that repetition gets large relative to knn, which
      // is why the fixture keeps the donors small (72 cells each) rather than raising knn.

      // Pseudotime dynamics
      obs_pseudotime: "palantir_pseudotime",
      n_splines:      4,
      dynamics_lam:   0.6,
      uns_dynamics:   "dynamics",

      // Cellular communities (small n for test speed)
      n_communities:       2,
      communities_alpha:   0.5,
      communities_method:  "hierarchical",

      // Trait associations. amyloid and braak follow the simulated severity, age does not,
      // so the test can check that the true associations are found and the null one is not.
      traits_csv:            resources_test.resolve("beyond_test_data/traits.csv"),
      trait_columns:         ["amyloid", "braak", "age"],
      proportion_transform:  "clr",

      // Pathway enrichment (local GMT file so no internet required in CI)
      de_results_csv:          resources_test.resolve("beyond_test_data/de_ExN.csv"),
      // gene names are the row index of the CSV (no named gene_column needed)
      gene_sets_file:          [resources_test.resolve("beyond_test_data/gene_sets.gmt")],
      pathway_method:          "prerank",
      output_pathway_csv:      "pathway_enrichment.csv",
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

    // the association table leaves the workflow as a CSV, not inside the h5mu
    assert state.containsKey("output_trait_associations_csv") :
      "State should contain key 'output_trait_associations_csv'. Found: ${state.keySet()}"
    assert state.output_trait_associations_csv.isFile() :
      "'output_trait_associations_csv' should be a file."
    def assoc_lines = state.output_trait_associations_csv.readLines()
    assert assoc_lines[0] ==
      "response,predictor,term,beta,se,stat,p_value,fdr_q,n,model,converged" :
      "Unexpected association CSV header: ${assoc_lines[0]}"
    // 3 traits x the subpopulations of the test atlas
    assert assoc_lines.size() > 2 :
      "Association CSV should hold more than one result row. Found: ${assoc_lines.size() - 1}"

    // The fixture plants a severity gradient that amyloid and braak follow and age does
    // not; if a step in the chain stops working, these flip.
    def assoc_header = assoc_lines[0].split(",")
    def i_predictor  = assoc_header.findIndexOf { it == "predictor" }
    def i_fdr        = assoc_header.findIndexOf { it == "fdr_q" }
    def min_q = [:].withDefault { 1.0 }
    assoc_lines.drop(1).each { line ->
      def f = line.split(",")
      def predictor = f[i_predictor]
      def q = f[i_fdr] as Double
      if (q < min_q[predictor]) { min_q[predictor] = q }
    }
    assert min_q["amyloid"] < 0.05 :
      "amyloid follows the simulated trajectory but was not detected (min fdr_q ${min_q['amyloid']})"
    assert min_q["braak"] < 0.05 :
      "braak follows the simulated trajectory but was not detected (min fdr_q ${min_q['braak']})"
    assert min_q["age"] > 0.05 :
      "age is independent of the simulated trajectory but came out significant (min fdr_q ${min_q['age']})"

    // enrichment results are one table too, not a directory of files
    assert state.containsKey("output_pathway_csv") :
      "State should contain key 'output_pathway_csv'. Found: ${state.keySet()}"
    assert state.output_pathway_csv.isFile() :
      "'output_pathway_csv' should be a file."
    def pathway_lines = state.output_pathway_csv.readLines()
    assert pathway_lines[0].startsWith("gene_set_library,method,") :
      "Unexpected enrichment CSV header: ${pathway_lines[0]}"
    assert pathway_lines.size() > 1 :
      "Enrichment CSV should hold at least one result row."

    // The DE tables are computed from the simulated counts, so the planted DISEASE_UP set
    // must come out enriched with a positive score.
    def pw_header = pathway_lines[0].split(",")
    def i_term    = pw_header.findIndexOf { it == "Term" }
    def i_nes     = pw_header.findIndexOf { it == "NES" }
    def disease_up = pathway_lines.drop(1).find { it.split(",")[i_term] == "DISEASE_UP" }
    assert disease_up != null :
      "DISEASE_UP missing from the enrichment results"
    assert (disease_up.split(",")[i_nes] as Double) > 0 :
      "DISEASE_UP should be positively enriched in the simulated DE results"

    "Output: $output"
  }
}
