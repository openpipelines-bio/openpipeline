nextflow.enable.dsl=2

include { trajectory_analysis } from params.rootDir + "/target/nextflow/workflows/beyond/trajectory_analysis/main.nf"
include { trajectory_analysis_test } from params.rootDir + "/target/_test/nextflow/test_workflows/beyond/trajectory_analysis_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

// Two fixtures, because neither answers the other's question.
//
//   beyond_test_data            PsychAD RADC, 152 real donors. Says whether the workflow
//                               survives a real annotation hierarchy, real sparsity and
//                               real donor metadata. It carries no compositional signal
//                               that survives multiple testing, so it can only be
//                               asserted on structurally.
//   beyond_simulated_test_data  12 donors with a planted severity gradient. Says whether
//                               the workflow *recovers* what is in the data, which is
//                               only checkable when we put it there.
workflow test_wf {
  resources_test = file(params.resources_test)

  Channel.fromList([
    [
      id:    "beyond_radc",
      input: resources_test.resolve("beyond_test_data/atlas.h5mu"),

      output:                         "beyond_radc.h5mu",
      output_trait_associations_csv:  "trait_associations.csv",

      // Atlas metadata. Subpopulation prevalence is taken within cell class, as the
      // BEYOND reference implementation does.
      obs_group:              "participant_id",
      obs_label:              "subpopulation",
      obs_normalize_within:   "cell_class",

      // Cellular landscape, at the reference parameters for the full cohort
      // (1.create.cellular.landscape.R: n_components=3, k=10, a=40).
      phate_n_components: 3,
      phate_knn:          10,
      phate_decay:        40,

      // Palantir: the root is picked from the donors without an AD diagnosis. The two
      // endpoints are pinned because Palantir's automatic terminal-state detection finds
      // none on this landscape - there is no branching structure in it to find.
      start_cluster_column:     "AD_status",
      start_cluster:            "No",
      terminal_states:          ["Donor_1007", "Donor_1024"],
      num_waypoints:            50,
      palantir_waypoint_knn:    15,
      palantir_n_components:    5,
      palantir_knn:             15,
      pseudotime_column:        "palantir_pseudotime",

      dynamics_lam:   0.6,

      n_communities:                    6,
      communities_alpha:                0.5,
      communities_correlation_method:   "spearman",
      communities_method:               "hierarchical",

      traits_csv:            resources_test.resolve("beyond_test_data/traits.csv"),
      trait_columns:         ["AD_status", "sex", "age"],
      proportion_transform:  "sqrt",
    ],
    [
      id:    "beyond_simulated",
      input: resources_test.resolve("beyond_simulated_test_data/atlas.h5mu"),

      output:                         "beyond_simulated.h5mu",
      output_trait_associations_csv:  "trait_associations.csv",

      obs_group: "participant_id",
      obs_label: "subpopulation",

      // 12 donors, so every neighbourhood parameter has to stay below that.
      phate_n_components: 3,
      phate_knn:          5,

      // Root is the healthy end of the simulated severity gradient, endpoints the severe
      // end. Pinned for the same reason as above, plus Palantir's automatic detection
      // solves a 10-eigenvector problem that 12 donors cannot support.
      start_id:                 "donor_01",
      terminal_states:          ["donor_11", "donor_12"],
      num_waypoints:            12,
      palantir_waypoint_knn:    5,
      palantir_n_components:    3,
      palantir_knn:             5,
      pseudotime_column:        "palantir_pseudotime",

      dynamics_lam:   0.6,

      n_communities:                    2,
      communities_alpha:                0.5,
      communities_correlation_method:   "spearman",
      communities_method:               "hierarchical",

      // amyloid and braak follow the simulated severity, age does not, so the test can
      // check that the true associations are found and the null one is not.
      traits_csv:            resources_test.resolve("beyond_simulated_test_data/traits.csv"),
      trait_columns:         ["amyloid", "braak", "age"],
      proportion_transform:  "clr",
    ]
  ])
  | map { state -> [state.id, state] }
  | trajectory_analysis
  | view { output ->
    assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

    def id = output[0]
    def state = output[1]
    assert state instanceof Map : "State should be a map. Found: ${state}"

    // -- shared structure ------------------------------------------------------------
    // The donor-level h5mu itself is checked by trajectory_analysis_test below.
    assert state.containsKey("output") : "Output should contain key 'output'."
    assert state.output.isFile() : "'output' should be a file."
    assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

    // The association table leaves the workflow as a CSV, not inside the h5mu
    assert state.containsKey("output_trait_associations_csv") :
      "State should contain key 'output_trait_associations_csv'. Found: ${state.keySet()}"
    assert state.output_trait_associations_csv.isFile() :
      "'output_trait_associations_csv' should be a file."
    def assoc_lines = state.output_trait_associations_csv.readLines()
    assert assoc_lines[0] ==
      "response,predictor,term,beta,se,stat,p_value,fdr_q,n,model,converged" :
      "Unexpected association CSV header: ${assoc_lines[0]}"
    assert assoc_lines.size() > 2 :
      "Association CSV should hold more than one result row. Found: ${assoc_lines.size() - 1}"

    // -- what only the simulated fixture can check -------------------------------------
    //
    // The real cohort carries no compositional association that survives BH correction
    // (152 donors; closest is AD_status x OPC at class level, q = 0.053), so only the
    // planted gradient can tell a working chain from a broken one that still emits
    // correctly shaped tables.
    if (id == "beyond_simulated") {
      def assoc_header = assoc_lines[0].split(",")
      def i_predictor  = assoc_header.findIndexOf { it == "predictor" }
      def i_fdr        = assoc_header.findIndexOf { it == "fdr_q" }
      def min_q = [:].withDefault { 1.0 }
      assoc_lines.drop(1).each { line ->
        def f = line.split(",")
        def q = f[i_fdr] as Double
        if (q < min_q[f[i_predictor]]) { min_q[f[i_predictor]] = q }
      }
      assert min_q["amyloid"] < 0.05 :
        "amyloid follows the simulated trajectory but was not detected (min fdr_q ${min_q['amyloid']})"
      assert min_q["braak"] < 0.05 :
        "braak follows the simulated trajectory but was not detected (min fdr_q ${min_q['braak']})"
      assert min_q["age"] > 0.05 :
        "age is independent of the simulated trajectory but came out significant (min fdr_q ${min_q['age']})"
    }

    "Output: $output"
  }
  | trajectory_analysis_test.run(
      fromState: { id, state ->
        def radc = id == "beyond_radc"
        [
          "input":                       state.output,
          "expected_n_obs":              radc ? 152 : 12,
          "expected_n_vars":             radc ? 65 : 9,
          "expected_n_communities":      radc ? 6 : 2,
          "expected_n_terminal_states":  2,
          "expected_n_components":       3,
          // Within-class normalisation is only used on the real fixture.
          "obs_normalize_within":        radc ? "cell_class" : null,
        ]
      }
    )
}
