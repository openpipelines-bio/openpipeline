workflow run_wf {
  take:
    input_ch

  main:

    // QC metrics are only calculated for the modalities listed here. Each entry
    // maps a modality to:
    //   - layer: the state key holding that modality's input layer.
    //   - var_qc_metrics_keys: state keys pointing to boolean .var columns
    //     (detected during single-sample processing) for which per-cell
    //     proportions are calculated. Empty for modalities without such columns.
    def qc_metrics_modalities = [
      "rna": [
        "layer": "rna_layer",
        "var_qc_metrics_keys": ["var_name_mitochondrial_genes", "var_name_ribosomal_genes"],
      ],
      "prot": [
        "layer": "prot_layer",
        "var_qc_metrics_keys": [],
      ],
    ].asImmutable()

    output_ch = input_ch
      // Make sure there is not conflict between the output from this workflow
      // And the output from any of the components
      | map {id, state ->
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
      }
      | process_singlesample_base.run(
        fromState: {id, state -> state},
        toState: [
          "input": "output",
          "modality": "output_modality"
        ]
      )
      // Calculate basic QC metrics for the relevant modalities (see
      // qc_metrics_modalities above). Mitochondrial and ribosomal genes were
      // already detected during single-sample processing, so they can be used
      // directly here.
      // Runs per modality before the modalities are merged back together.
      // process_singlesample_base assigns the same id to every modality of a
      // sample, so make the id unique per modality to let the qc component run
      // per modality. The original id is kept in _meta.join_id to restore it
      // before grouping the modalities back together.
      | map {id, state -> ["${id}_${state.modality}", state + ["_meta": ["join_id": id]]]}
      | calculate_qc_metrics.run(
        runIf: {id, state -> !state.skip_qc_metrics && qc_metrics_modalities.containsKey(state.modality)},
        fromState: {id, state ->
          def modality_config = qc_metrics_modalities[state.modality]
          // For modalities with boolean .var columns, a user-provided
          // --var_qc_metrics takes precedence; otherwise the configured columns
          // (e.g. mitochondrial/ribosomal for rna) are used when present.
          def var_qc_metrics = null
          if (modality_config.var_qc_metrics_keys) {
            def detected_var_qc_metrics = modality_config.var_qc_metrics_keys.collect{state[it]}.findAll{it != null}
            var_qc_metrics = state.var_qc_metrics ?: detected_var_qc_metrics
          }
          def newState = [
            "input": state.input,
            "modality": state.modality,
            "layer": state[modality_config.layer],
            "top_n_vars": state.top_n_vars,
            "log1p_transform": state.log1p_transform,
            "output_obs_num_nonzero_vars": state.output_obs_num_nonzero_vars,
            "output_obs_total_counts_vars": state.output_obs_total_counts_vars,
            "output_var_num_nonzero_obs": state.output_var_num_nonzero_obs,
            "output_var_total_counts_obs": state.output_var_total_counts_obs,
            "output_var_obs_mean": state.output_var_obs_mean,
            "output_var_pct_dropout": state.output_var_pct_dropout,
            "output_compression": "gzip"
          ]
          if (var_qc_metrics) {
            newState += ["var_qc_metrics": var_qc_metrics]
          }
          return newState
        },
        toState: ["input": "output"]
      )
      // Restore the original sample id so the modalities group back together.
      | map {id, state -> [state._meta.join_id, state]}
      | groupTuple(by: 0, sort: "hash")
      | map { id, states ->
          def new_input = states.collect{it.input}
          def modalities = states.collect{it.modality}.unique()
          def other_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
          }.minus(["output", "input", "modality"])
          def new_state = other_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "Arguments should be the same across modalities. Please report this \
                                                 as a bug. Argument name: $argument_name, \
                                                 argument value: $argument_values"
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
          }
          [id, new_state + ["input": new_input, "modalities": modalities]]
      }
      | view {"Input merge channel: $it"}
      | merge.run(
        fromState: [
          "input": "input",
          "output": "workflow_output"
        ],
        toState: ["output": "output"]
      )
      | view {"After smerging: $it"}
      | intersect_obs.run(
        runIf: {id, state -> state.intersect_obs},
        fromState:[
          "input": "output",
          "modalities": "modalities",
          "output": "workflow_output"
        ],
        toState: {id, output, state -> 
          ["output": output.output]
        }
      )
      | view {"After singlesample processing: $it"}
      | setState(["output"])
    // output_ch = modalities_ch

  emit:
    output_ch
}