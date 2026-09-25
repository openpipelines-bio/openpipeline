nextflow.enable.dsl=2

include { trajectory_analysis } from params.rootDir + "/target/nextflow/workflows/beyond/trajectory_analysis/main.nf"

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

      output:                     "beyond_radc.h5mu",
      output_proportions_csv:     "proportions.csv",
      output_phate_csv:           "phate.csv",
      output_pseudotime_csv:      "pseudotime.csv",
      output_dynamics_csv:        "dynamics.csv",
      output_dynamics_stats_csv:  "dynamics_stats.csv",
      output_communities_csv:     "communities.csv",

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

      output:                     "beyond_simulated.h5mu",
      output_proportions_csv:     "proportions.csv",
      output_phate_csv:           "phate.csv",
      output_pseudotime_csv:      "pseudotime.csv",
      output_dynamics_csv:        "dynamics.csv",
      output_dynamics_stats_csv:  "dynamics_stats.csv",
      output_communities_csv:     "communities.csv",

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
    assert state.containsKey("output") : "Output should contain key 'output'."
    assert state.output.isFile() : "'output' should be a file."
    assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

    // Steps 2-5 are table in / table out: the landscape never enters the h5mu.
    [
      "output_proportions_csv",
      "output_phate_csv",
      "output_pseudotime_csv",
      "output_dynamics_csv",
      "output_dynamics_stats_csv",
      "output_communities_csv",
    ].each { key ->
      assert state.containsKey(key) : "State should contain key '${key}'. Found: ${state.keySet()}"
      assert state[key].isFile() : "'${key}' should be a file."
    }

    // One row per donor in the landscape tables - not one per cell, which is the whole
    // point of the rewrite: 152 donors out of 32377 nuclei, 12 out of 864.
    def n_donors = id == "beyond_radc" ? 152 : 12
    def n_labels = id == "beyond_radc" ? 65 : 9
    def n_communities = id == "beyond_radc" ? 6 : 2

    def phate_lines = state.output_phate_csv.readLines()
    assert phate_lines[0] == "participant_id,phate_1,phate_2,phate_3" :
      "Unexpected PHATE CSV header: ${phate_lines[0]}"
    assert phate_lines.size() == n_donors + 1 :
      "PHATE table should hold ${n_donors} donors, found ${phate_lines.size() - 1}"

    def pt_lines = state.output_pseudotime_csv.readLines()
    assert pt_lines[0].startsWith("participant_id,palantir_pseudotime,palantir_entropy") :
      "Unexpected pseudotime CSV header: ${pt_lines[0]}"
    assert pt_lines.size() == n_donors + 1 :
      "Pseudotime table should hold ${n_donors} donors, found ${pt_lines.size() - 1}"
    // Terminal states were pinned, so there must be one fate column per endpoint
    def n_fate = pt_lines[0].split(",").count { it.startsWith("fate_") }
    assert n_fate == 2 :
      "Expected 2 fate columns for the 2 pinned terminal states, found ${n_fate}"

    def comm_lines = state.output_communities_csv.readLines()
    assert comm_lines[0] == "label,community_id" :
      "Unexpected communities CSV header: ${comm_lines[0]}"
    assert comm_lines.size() == n_labels + 1 :
      "Communities table should hold ${n_labels} subpopulations, found ${comm_lines.size() - 1}"
    def comm_ids = comm_lines.drop(1).collect { it.split(",")[1] }.toSet()
    assert comm_ids.size() == n_communities :
      "Expected ${n_communities} communities, found ${comm_ids.size()}: ${comm_ids}"

    def dyn_lines = state.output_dynamics_csv.readLines()
    assert dyn_lines[0] == "label,pseudotime,proportion_fitted" :
      "Unexpected dynamics CSV header: ${dyn_lines[0]}"

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

    // -- what only the real fixture can check ------------------------------------------
    //
    // Within-class normalisation: a donor's row sums to the number of cell classes they
    // have nuclei in (8 here), not to 1. A global normalisation would give 1.
    if (id == "beyond_radc") {
      def prop_lines = state.output_proportions_csv.readLines()
      assert prop_lines.size() == n_donors + 1 :
        "Proportion table should hold ${n_donors} donors, found ${prop_lines.size() - 1}"
      assert prop_lines[0].split(",").size() == n_labels + 1 :
        "Proportion table should hold ${n_labels} subpopulations, found ${prop_lines[0].split(",").size() - 1}"
      def row_sum = prop_lines[1].split(",").drop(1).collect { it as Double }.sum()
      assert row_sum > 7.0 && row_sum < 8.01 :
        "Within-class proportions should sum to the number of classes present (<= 8), found ${row_sum}"
    }

    "Output: $output"
  }
}
