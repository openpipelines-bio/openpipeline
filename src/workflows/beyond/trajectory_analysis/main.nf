// beyond/trajectory_analysis - atlas h5mu -> full BEYOND outputs
//
// Input:  single channel event [ id, state ] where state.input is the atlas h5mu
//         carrying obs[--obs_group] + obs[--obs_label] (see beyond/subpopulation_clustering).
// Output: the atlas h5mu with the proportion matrix in uns, plus the group-level
//         result tables as CSV files.
//
// Step 1 is the only step that reads cell-level data. It converts the atlas into a
// group x label proportion table, and every step after it is table in / table out -
// the BEYOND cellular landscape is a donor landscape, not a cell landscape.
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

    // -- 1. Group x label proportion matrix ----------------------------------------
    //
    // The only step that opens the h5mu. The CSV it writes is what every downstream
    // step consumes.
    | calculate_label_proportions.run(
        fromState: { id, state -> [
          "input":       state.input,
          "obs_group":   state.obs_group,
          "obs_label":   state.obs_label,
          "obs_normalize_within": state.obs_normalize_within,
          "uns_output":  state.uns_proportions,
          "output":      state.workflow_output,
          "output_csv":  state.output_proportions_csv,
        ]},
        toState: { id, output, state ->
          state + [
            "input":                    output.output,
            "output_proportions_csv":   output.output_csv,
          ]
        }
      )

    // -- 2. PHATE cellular landscape (donor x subpopulation proportions) -----------
    | phate.run(
        fromState: { id, state -> [
          "input_table":  state.output_proportions_csv,
          "id_column":    state.obs_group,
          "output_table": state.output_phate_csv,
          "n_components": state.phate_n_components,
          "knn":          state.phate_knn,
          "decay":        state.phate_decay,
        ]},
        toState: { id, output, state ->
          state + [ "output_phate_csv": output.output_table ]
        }
      )

    // -- 3. Palantir pseudotime + fate probabilities on the landscape --------------
    //
    // The trait table doubles as the label source for --start_cluster, so a root
    // can be named as e.g. "the non-demented donors" rather than a single identifier.
    | palantir.run(
        fromState: { id, state -> [
          "input_table":             state.output_phate_csv,
          "id_column":               state.obs_group,
          "metadata":                state.traits_csv,
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
          "output_table":            state.output_pseudotime_csv,
        ]},
        toState: { id, output, state ->
          state + [ "output_pseudotime_csv": output.output_table ]
        }
      )

    // -- 4. Spline-fitted proportion dynamics along pseudotime ---------------------
    | fit_proportion_dynamics.run(
        fromState: { id, state -> [
          "input":              state.output_proportions_csv,
          "pseudotime":         state.output_pseudotime_csv,
          "id_column":          state.obs_group,
          "pseudotime_column":  state.pseudotime_column,
          "lam":                state.dynamics_lam,
          "output":             state.output_dynamics_csv,
          "output_stats":       state.output_dynamics_stats_csv,
        ]},
        toState: { id, output, state ->
          state + [
            "output_dynamics_csv":        output.output,
            "output_dynamics_stats_csv":  output.output_stats,
          ]
        }
      )

    // -- 5. Label communities (co-occurrence + dynamics) ----------------------------
    | label_communities.run(
        fromState: { id, state -> [
          "input":               state.output_proportions_csv,
          "dynamics":            state.output_dynamics_csv,
          "id_column":           state.obs_group,
          "n_communities":       state.n_communities,
          "alpha":               state.communities_alpha,
          "correlation_method":  state.communities_correlation_method,
          "method":              state.communities_method,
          "output":              state.output_communities_csv,
        ]},
        toState: { id, output, state ->
          state + [ "output_communities_csv": output.output ]
        }
      )

    // -- 6. Trait associations on the proportion table ------------------------------
    | test_associations.run(
        fromState: { id, state -> [
          "input":                state.output_proportions_csv,
          "metadata":             state.traits_csv,
          "input_join_column":    state.obs_group,
          "metadata_join_column": state.obs_group,
          "predictor_columns":    state.trait_columns,
          "covariate_columns":    state.covariate_columns,
          "random_effect_column": state.random_effect_column,
          "transform":            state.proportion_transform,
          "fdr_method":           state.fdr_method,
          "fdr_scope":            state.fdr_scope,
          "output":               state.output_trait_associations_csv,
        ]},
        toState: { id, output, state ->
          state + [ "output_trait_associations_csv": output.output ]
        }
      )

    | map { id, state ->
        def out = [
          "output":                         state.input,
          "output_proportions_csv":         state.output_proportions_csv,
          "output_phate_csv":               state.output_phate_csv,
          "output_pseudotime_csv":          state.output_pseudotime_csv,
          "output_dynamics_csv":            state.output_dynamics_csv,
          "output_dynamics_stats_csv":      state.output_dynamics_stats_csv,
          "output_communities_csv":         state.output_communities_csv,
          "output_trait_associations_csv":  state.output_trait_associations_csv,
        ]
        [ id, out ]
      }

  emit:
  output_ch
}
