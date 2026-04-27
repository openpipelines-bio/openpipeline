// beyond/trajectory_analysis - atlas h5mu -> full BEYOND outputs
//
// Input:  single channel event [ id, state ] where state.input is the atlas h5mu
//         produced by beyond/atlas_building (obs["participant_id"] + obs["subpopulation"]).
// Output: h5mu with all BEYOND annotations + optional CSV files.
//
// Channel convention throughout: [ id, state_map ]

workflow run_wf {
  take:
  input_ch  // [ id, { input: atlas.h5mu, ..params.. } ]

  main:

  // Preserve the final output filename so it is available after each toState swap.
  prep_ch = input_ch
    | map { id, state ->
        [ id, state + [ "workflow_output": state.output ] ]
      }

  // -- 1. Participant x subpopulation proportion matrix --------------------------
  prop_ch = prep_ch
    | calculate_proportions.run(
        fromState: { id, state -> [
          "input":              state.input,
          "obs_participant_id": state.obs_participant_id,
          "obs_subpopulation":  state.obs_subpopulation,
          "uns_output":         state.uns_proportions,
          "obsm_output":        state.uns_proportions,
          "output":             state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 2. PHATE cellular landscape (input: proportion matrix in obsm) -------------
  phate_ch = prop_ch
    | phate.run(
        fromState: { id, state -> [
          "input":      state.input,
          "obsm_input": state.uns_proportions,
          "output":     state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 3. Palantir pseudotime + fate probabilities (input: X_phate) ---------------
  palantir_ch = phate_ch
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

  // -- 4. VIA pseudotime (alternative trajectory, runs on same X_phate) ----------
  via_ch = palantir_ch
    | via.run(
        fromState: { id, state -> [
          "input":       state.input,
          "obsm_key":    "X_phate",
          "obs_cluster": state.obs_cluster,
          "root_user":   state.root_user,
          "knn":         state.via_knn,
          "output":      state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 5. GAM-fitted proportion dynamics along pseudotime ------------------------
  dynamics_ch = via_ch
    | pseudotime_dynamics.run(
        fromState: { id, state -> [
          "input":                state.input,
          "obs_pseudotime":       state.obs_pseudotime,
          "obs_participant_id":   state.obs_participant_id,
          "uns_proportions":      state.uns_proportions,
          "n_splines":            state.n_splines,
          "lam":                  state.dynamics_lam,
          "uns_output":           state.uns_dynamics,
          "output":               state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 6. Cellular community detection (co-occurrence + dynamics) ----------------
  communities_ch = dynamics_ch
    | cellular_communities.run(
        fromState: { id, state -> [
          "input":             state.input,
          "obs_subpopulation": state.obs_subpopulation,
          "uns_proportions":   state.uns_proportions,
          "uns_dynamics":      state.uns_dynamics,
          "n_communities":     state.n_communities,
          "alpha":             state.communities_alpha,
          "method":            state.communities_method,
          "output":            state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

  // -- 7. Linear mixed-model trait associations -----------------------------------
  traits_ch = communities_ch
    | trait_associations.run(
        fromState: { id, state -> [
          "input":                 state.input,
          "uns_proportions":       state.uns_proportions,
          "traits_csv":            state.traits_csv,
          "participant_id_column": state.obs_participant_id,
          "trait_columns":         state.trait_columns,
          "covariate_columns":     state.covariate_columns,
          "random_effect_column":  state.random_effect_column,
          "output":                state.workflow_output,
          "output_csv":            state.output_trait_associations_csv,
        ]},
        toState: { id, output, state ->
          state + [
            "input":                          output.output,
            "output_trait_associations_csv":  output.output_csv,
          ]
        }
      )

  // -- 8. Pathway enrichment (optional - skipped when de_results_csv is null) -----
  //
  // Branch into two sub-channels: events with DE results run pathway_enrichment;
  // events without skip directly to the final setState.

  (with_de_ch, no_de_ch) = traits_ch
    | branch { id, state ->
        with_de: state.de_results_csv != null
        no_de:   true
      }

  enriched_ch = with_de_ch
    | pathway_enrichment.run(
        fromState: { id, state -> [
          "input":          state.input,
          "input_degenes":  state.de_results_csv,
          "gene_column":    state.gene_column,
          "gene_sets":      state.gene_sets,
          "method":         state.pathway_method,
          "output":         state.workflow_output,
          "output_csv_dir": state.output_pathway_csv_dir ?: "pathway_enrichment_results",
        ]},
        toState: { id, output, state ->
          state + [
            "input":                    output.output,
            "output_pathway_csv_dir":   output.output_csv_dir,
          ]
        }
      )

  output_ch = enriched_ch
    .mix(no_de_ch)
    | map { id, state -> [ id, [ "output": state.input ] ] }

  emit:
  output_ch
}
