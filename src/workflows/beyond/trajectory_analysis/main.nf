// beyond/trajectory_analysis - atlas h5mu -> full BEYOND outputs
//
// Input:  single channel event [ id, state ] where state.input is the atlas h5mu
//         carrying obs[--obs_group] + obs[--obs_label] (see beyond/subpopulation_clustering).
// Output: h5mu with the BEYOND annotations, plus the association table and the optional
//         enrichment table as CSV files.
//
// Steps 1-5 pass the h5mu along; steps 6-7 work on tables and leave the h5mu untouched.
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
    // The CSV copy is what the (MuData-free) association step downstream consumes.
    | calculate_label_proportions.run(
        fromState: { id, state -> [
          "input":       state.input,
          "obs_group":   state.obs_group,
          "obs_label":   state.obs_label,
          "uns_output":  state.uns_proportions,
          "obsm_output": state.uns_proportions,
          "output":      state.workflow_output,
          "output_csv":  "${id}.proportions.csv",
        ]},
        toState: { id, output, state ->
          state + [
            "input":            output.output,
            "proportions_csv":  output.output_csv,
          ]
        }
      )

    // -- 2. PHATE cellular landscape (input: proportion matrix in obsm) -------------
    | phate.run(
        fromState: { id, state -> [
          "input":      state.input,
          "obsm_input": state.uns_proportions,
          "output":     state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

    // -- 3. Palantir pseudotime + fate probabilities (input: X_phate) ---------------
    | palantir.run(
        fromState: { id, state -> [
          "input":                   state.input,
          "obsm_input":              "X_phate",
          "start_cell":              state.start_cell,
          "start_cell_cluster":      state.start_cell_cluster,
          "start_cell_obs_key":      state.start_cell_obs_key,
          "terminal_states_obs_key": state.terminal_states_obs_key,
          "num_waypoints":           state.num_waypoints,
          "n_components":            state.palantir_n_components,
          "knn":                     state.palantir_knn,
          "output":                  state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

    // -- 4. Spline-fitted proportion dynamics along pseudotime ---------------------
    | fit_proportion_dynamics.run(
        fromState: { id, state -> [
          "input":           state.input,
          "obs_pseudotime":  state.obs_pseudotime,
          "obs_group":       state.obs_group,
          "uns_proportions": state.uns_proportions,
          "n_splines":       state.n_splines,
          "lam":             state.dynamics_lam,
          "uns_output":      state.uns_dynamics,
          "output":          state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

    // -- 5. Cellular community detection (co-occurrence + dynamics) ----------------
    | label_communities.run(
        fromState: { id, state -> [
          "input":           state.input,
          "obs_label":       state.obs_label,
          "uns_proportions": state.uns_proportions,
          "uns_dynamics":    state.uns_dynamics,
          "n_communities":   state.n_communities,
          "alpha":           state.communities_alpha,
          "method":          state.communities_method,
          "output":          state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

    // -- 6. Trait associations on the proportion table ------------------------------
    //
    // Tabular in / tabular out: the proportion CSV from step 1 is joined with the trait
    // table on the group column, so the h5mu (state.input) is left untouched here.
    | test_associations.run(
        fromState: { id, state -> [
          "input":                state.proportions_csv,
          "metadata":             state.traits_csv,
          "join_on":              state.obs_group,
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

    // -- 7. Pathway enrichment (skipped when no DE results are given) ---------------
    //
    // Also CSV in / CSV out, so it needs no branch/mix: runIf skips the step and leaves
    // the state as it was.
    | gseapy.run(
        runIf: { id, state -> state.de_results_csv != null },
        fromState: { id, state -> [
          "input":          state.de_results_csv,
          "gene_column":    state.gene_column,
          "gene_sets":      state.gene_sets,
          "gene_sets_file": state.gene_sets_file,
          "method":         state.pathway_method,
          "output":         state.output_pathway_csv,
        ]},
        toState: { id, output, state ->
          state + [ "output_pathway_csv": output.output ]
        }
      )

    | map { id, state ->
        // The association table always exists; the enrichment table only when gseapy ran.
        def out = [
          "output":                         state.input,
          "output_trait_associations_csv":  state.output_trait_associations_csv,
        ]
        if (state.de_results_csv != null) {
          out = out + [ "output_pathway_csv": state.output_pathway_csv ]
        }
        [ id, out ]
      }

  emit:
  output_ch
}
