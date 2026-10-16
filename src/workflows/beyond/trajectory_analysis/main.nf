// beyond/trajectory_analysis - atlas h5mu -> donor-level h5mu with the BEYOND results
//
// Input:  single channel event [ id, state ] where state.input is the atlas h5mu
//         carrying obs[--obs_group] + obs[--obs_label] (see beyond/subpopulation_clustering).
// Output: a donor-level h5mu (one observation per donor, one variable per
//         subpopulation) plus the trait associations as a CSV.
//
// Step 1 is the only step that reads cell-level data. It turns the atlas into a new,
// donor-level h5mu, and every step after it reads and extends that object - the
// BEYOND cellular landscape is a donor landscape, not a cell landscape.
//
// Channel convention throughout: [ id, state_map ]

workflow run_wf {
  take:
  input_ch  // [ id, { input: atlas.h5mu, ..params.. } ]

  main:

  output_ch = input_ch

    // Preserve the final output filename so it is available after each toState swap.
    | map { id, state ->
        [ id, state + [ "workflow_output": state.output ] ]
      }

    // -- 1. Donor x subpopulation proportions, as a new donor-level h5mu -----------
    | calculate_obs_label_proportions.run(
        fromState: { id, state -> [
          "input":                state.input,
          "obs_group":            state.obs_group,
          "obs_label":            state.obs_label,
          "obs_normalize_within": state.obs_normalize_within,
          "output":               state.workflow_output,
        ]},
        args: [
          "output_modality": "proportions",
          "obsm_output":     "X_proportions",
        ],
        toState: [ "input": "output" ]
      )

    // -- 1b. Donor traits onto the donors ------------------------------------------
    | join_csv.run(
        fromState: { id, state -> [
          "input":     state.input,
          "input_csv": state.traits_csv,
          "obs_key":   state.obs_group,
          "csv_key":   state.obs_group,
          "output":    state.workflow_output,
        ]},
        args: [ "modality": "proportions" ],
        toState: [ "input": "output" ]
      )

    // -- 2. PHATE cellular landscape (donor x subpopulation proportions) -----------
    | phate.run(
        fromState: { id, state -> [
          "input":        state.input,
          "n_components": state.phate_n_components,
          "knn":          state.phate_knn,
          "decay":        state.phate_decay,
          "output":       state.workflow_output,
        ]},
        args: [
          "modality":    "proportions",
          "obsm_input":  "X_proportions",
          "obsm_output": "X_phate",
        ],
        toState: [ "input": "output" ]
      )

    // -- 3. Palantir pseudotime + fate probabilities on the landscape --------------
    //
    // The traits are in .obs, so a root can be named as e.g. "the non-demented
    // donors" (--start_cluster) rather than a single identifier.
    | palantir.run(
        fromState: { id, state -> [
          "input":                   state.input,
          "start_id":                state.start_id,
          "start_cluster":           state.start_cluster,
          "start_cluster_column":    state.start_cluster_column,
          "terminal_states":         state.terminal_states,
          "terminal_states_column":  state.terminal_states_column,
          "num_waypoints":           state.num_waypoints,
          "n_components":            state.palantir_n_components,
          "knn":                     state.palantir_knn,
          "waypoint_knn":            state.palantir_waypoint_knn,
          "pseudotime_column":       state.pseudotime_column,
          "output":                  state.workflow_output,
        ]},
        args: [
          "modality":   "proportions",
          "obsm_input": "X_phate",
        ],
        toState: [ "input": "output" ]
      )

    // -- 4. Spline-fitted proportion dynamics along pseudotime ---------------------
    | fit_proportion_dynamics.run(
        fromState: { id, state -> [
          "input":          state.input,
          "obs_pseudotime": state.pseudotime_column,
          "lam":            state.dynamics_lam,
          "output":         state.workflow_output,
        ]},
        args: [ "modality": "proportions" ],
        toState: [ "input": "output" ]
      )

    // -- 5. Subpopulation communities (co-occurrence + dynamics) --------------------
    | subpopulation_communities.run(
        fromState: { id, state -> [
          "input":               state.input,
          "n_communities":       state.n_communities,
          "alpha":               state.communities_alpha,
          "correlation_method":  state.communities_correlation_method,
          "method":              state.communities_method,
          "output":              state.workflow_output,
        ]},
        args: [ "modality": "proportions" ],
        toState: [ "input": "output" ]
      )

    // -- 6. Trait associations on the proportions ----------------------------------
    | test_associations.run(
        fromState: { id, state -> [
          "input":                state.input,
          "predictor_columns":    state.trait_columns,
          "covariate_columns":    state.covariate_columns,
          "random_effect_column": state.random_effect_column,
          "transform":            state.proportion_transform,
          "fdr_method":           state.fdr_method,
          "fdr_scope":            state.fdr_scope,
          "output":               state.output_trait_associations_csv,
        ]},
        args: [ "modality": "proportions" ],
        toState: [ "output_trait_associations_csv": "output" ]
      )

    | map { id, state ->
        [ id, [
          "output":                         state.input,
          "output_trait_associations_csv":  state.output_trait_associations_csv,
        ]]
      }

  emit:
  output_ch
}
