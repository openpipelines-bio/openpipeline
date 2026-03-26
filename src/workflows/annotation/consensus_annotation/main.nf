workflow run_wf {
  take:
    input_ch

  main:
    // Stash the final output path so sub-workflow runs don't claim it.
    // Also initialise current_target (the file progressively enriched by the
    // copy steps) to the original input; the first enabled copy step will
    // overwrite it.
    stashed_ch = input_ch
      | map { id, state ->
          assert !state.skip_scanvi || !state.skip_celltypist || !state.skip_singler :
              "At least one annotation method must be enabled. " +
              "Use --skip_scanvi / --skip_celltypist / --skip_singler to disable individual methods."
          [id, state + [
            "workflow_output":  state.output,
            "current_target":   state.input
          ]]
        }

    // -----------------------------------------------------------------------
    // Run enabled annotation methods in parallel.
    // runIf lets the event pass through unchanged when a method is skipped,
    // so all three channels always emit exactly one event per sample — this
    // keeps the downstream .join() safe regardless of which methods are on.
    // -----------------------------------------------------------------------

    scanvi_ch = stashed_ch
      | scanvi_scarches_workflow.run(
          runIf: { id, state -> !state.skip_scanvi },
          fromState: { id, state -> [
            "id":                          id,
            "input":                       state.input,
            "modality":                    state.modality,
            "layer":                       state.scanvi_layer,
            "input_obs_batch_label":       state.scanvi_input_obs_batch_label,
            "reference":                   state.scanvi_reference,
            "reference_obs_target":        state.scanvi_reference_obs_target,
            "reference_obs_batch_label":   state.scanvi_reference_obs_batch_label,
            "reference_var_hvg":           state.scanvi_reference_var_hvg,
            "output_obs_predictions":      state.scanvi_output_obs_predictions,
            "output_obs_probability":      state.scanvi_output_obs_probability,
            "output_obsm_integrated":      state.scanvi_output_obsm_integrated,
            "output":                      "scanvi_annotated.h5mu"
          ]},
          toState: { id, output, state ->
            state + ["scanvi_output": output.output]
          }
        )

    celltypist_ch = stashed_ch
      | celltypist_workflow.run(
          runIf: { id, state -> !state.skip_celltypist },
          fromState: { id, state -> [
            "input":                       state.input,
            "modality":                    state.modality,
            "input_layer":                 state.celltypist_layer,
            "reference":                   state.celltypist_reference,
            "reference_obs_target":        state.celltypist_reference_obs_target,
            "model":                       state.celltypist_model,
            "majority_voting":             state.celltypist_majority_voting,
            "output_obs_predictions":      state.celltypist_output_obs_predictions,
            "output_obs_probability":      state.celltypist_output_obs_probability,
            "output":                      "celltypist_annotated.h5mu"
          ]},
          toState: { id, output, state ->
            state + ["celltypist_output": output.output]
          }
        )

    singler_ch = stashed_ch
      | singler.run(
          runIf: { id, state -> !state.skip_singler },
          fromState: { id, state -> [
            "input":                           state.input,
            "modality":                        state.modality,
            "input_layer":                     state.singler_layer,
            "reference":                       state.singler_reference,
            "reference_obs_target":            state.singler_reference_obs_target,
            "output_obs_predictions":          state.singler_output_obs_predictions,
            "output_obs_probability":          state.singler_output_obs_probability,
            "output_obs_delta_next":           state.singler_output_obs_delta_next,
            "output_obs_pruned_predictions":   state.singler_output_obs_pruned_predictions,
            "output_obsm_scores":              state.singler_output_obsm_scores,
            "output":                          "singler_annotated.h5mu"
          ]},
          args: ["prune": true],
          toState: { id, output, state ->
            state + ["singler_output": output.output]
          }
        )

    // -----------------------------------------------------------------------
    // Join all three branches.  Because runIf passes events through when a
    // method is skipped, every sample always has exactly one event per
    // channel.  Skipped methods simply have no *_output key in their state.
    // -----------------------------------------------------------------------

    output_ch = scanvi_ch
      .join(celltypist_ch, failOnMismatch: true, failOnDuplicate: true)
      .join(singler_ch, failOnMismatch: true, failOnDuplicate: true)
      | map { id, scanvi_state, celltypist_state, singler_state ->
          // Use the scanvi branch as the base state (carries all shared keys
          // including workflow_output and current_target) and pull the output
          // paths from the other two branches.
          [id, scanvi_state + [
            "celltypist_output": celltypist_state.celltypist_output,
            "singler_output":    singler_state.singler_output
          ]]
        }

      // ---------------------------------------------------------------------
      // Copy each enabled method's annotation slots back onto the original
      // input file.  current_target starts as state.input (set above) and is
      // updated to the output of each copy step, so the three steps chain
      // sequentially without overwriting each other.
      // ---------------------------------------------------------------------

      | copy_h5mu_slots.run(
          runIf: { id, state -> !state.skip_scanvi },
          fromState: { id, state -> [
            "input":    state.scanvi_output,
            "target":   state.current_target,
            "modality": state.modality,
            "obs": [
              state.scanvi_output_obs_predictions,
              state.scanvi_output_obs_probability,
              "scanvi_integration_leiden"
            ],
            "obsm": [
              state.scanvi_output_obsm_integrated,
              "X_scanvi_umap"
            ],
            "uns":  ["scanvi_integration_neighbors"],
            "obsp": ["scanvi_integration_distances", "scanvi_integration_connectivities"],
            "output": "annotated.h5mu"
          ]},
          toState: { id, output, state ->
            state + ["current_target": output.output]
          }
        )

      | copy_h5mu_slots.run(
          runIf: { id, state -> !state.skip_celltypist },
          fromState: { id, state -> [
            "input":    state.celltypist_output,
            "target":   state.current_target,
            "modality": state.modality,
            "obs": [
              state.celltypist_output_obs_predictions,
              state.celltypist_output_obs_probability
            ],
            "output": "annotated.h5mu"
          ]},
          toState: { id, output, state ->
            state + ["current_target": output.output]
          }
        )

      | copy_h5mu_slots.run(
          runIf: { id, state -> !state.skip_singler },
          fromState: { id, state -> [
            "input":    state.singler_output,
            "target":   state.current_target,
            "modality": state.modality,
            "obs": [
              state.singler_output_obs_predictions,
              state.singler_output_obs_probability,
              state.singler_output_obs_delta_next,
              state.singler_output_obs_pruned_predictions
            ],
            "obsm": [state.singler_output_obsm_scores],
            "output": "annotated.h5mu"
          ]},
          toState: { id, output, state ->
            state + ["current_target": output.output]
          }
        )

      // ---------------------------------------------------------------------
      // Consensus vote.  Build prediction and probability lists dynamically
      // from the methods that were actually run.
      // ---------------------------------------------------------------------
      | consensus_vote.run(
          fromState: { id, state ->
            def preds = []
            def probs = []
            if (!state.skip_scanvi) {
              preds.add(state.scanvi_output_obs_predictions)
              probs.add(state.scanvi_output_obs_probability)
            }
            if (!state.skip_celltypist) {
              preds.add(state.celltypist_output_obs_predictions)
              probs.add(state.celltypist_output_obs_probability)
            }
            if (!state.skip_singler) {
              preds.add(state.singler_output_obs_predictions)
              probs.add(state.singler_output_obs_probability)
            }
            [
              "input":                  state.current_target,
              "modality":               state.modality,
              "input_obs_predictions":  preds,
              "input_obs_probabilities": probs,
              "weights":                state.consensus_weights,
              "tie_label":              state.consensus_tie_label,
              "output_obs_predictions": state.output_obs_predictions,
              "output_obs_score":       state.output_obs_score,
              "output":                 state.workflow_output
            ]
          },
          toState: ["output": "output"]
        )
      | setState(["output"])

  emit:
    output_ch
}
