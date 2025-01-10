workflow run_wf {
  take:
  input_ch

  main:
    preproc_ch = input_ch
    // Avoid conflict between output from component and output for this workflow
    | map {id, state -> 
      assert state.output, "Output must be defined"
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // Add default for var_qc_metrics component
    | map {id, state ->
      def var_qc_default = []
      // Remove the var_qc_metric argument from the state if its value is null (not specified)
      def new_state = state.findAll { it.key != "var_qc_metrics" || it.value == null }
      if (state.var_name_mitochondrial_genes) {
        var_qc_default.add(state.var_name_mitochondrial_genes)
      }
      // Get the new state, but make sure to overwrite var_qc_metrics if the user has set it.
      new_state = ["var_qc_metrics": var_qc_default.join(",")] + new_state
      [id, new_state]
    }

    with_grep_ch = preproc_ch
    | filter { it -> it[1].var_name_mitochondrial_genes }
    | grep_annotation_column.run(
      fromState: { id, state ->
        def stateMapping = [
          "input": state.input,
          "modality": state.modality,
          "input_column": state.var_gene_names,
          "matrix": "var",
          "output_match_column": state.var_name_mitochondrial_genes,
          "regex_pattern": state.mitochondrial_gene_regex,
          "input_layer": state.layer,
        ]
        stateMapping.output_fraction_column = state.obs_name_mitochondrial_fraction ? state.obs_name_mitochondrial_fraction: "fraction_$state.var_name_mitochondrial_genes"
        return stateMapping
      },
      toState: ["input": "output"]
    )

    without_grep_ch = preproc_ch
      | filter { it -> !it[1].var_name_mitochondrial_genes }

    output_ch = without_grep_ch.mix(with_grep_ch)
      | calculate_qc_metrics.run(
        fromState: { id, state ->
          def newState = [
            "input": state.input,
            "modality": state.modality,
            "layer": state.layer,
            // TODO: remove this workaround when Viash issue is resolved:
            //       'top_n_vars': list(map(int, r''.split(';'))),
            //     ValueError: invalid literal for int() with base 10: ''
            // See https://github.com/viash-io/viash/issues/619
            "top_n_vars": state.top_n_vars ? state.top_n_vars : null,
            "var_qc_metrics_fill_na_value": state.var_qc_metrics_fill_na_value,
            "output_obs_num_nonzero_vars": state.output_obs_num_nonzero_vars,
            "output_obs_total_counts_vars": state.output_obs_total_counts_vars,
            "output_var_num_nonzero_obs": state.output_var_num_nonzero_obs,
            "output_var_total_counts_obs": state.output_var_total_counts_obs,
            "output_var_obs_mean": state.output_var_obs_mean,
            "output_var_pct_dropout": state.output_var_pct_dropout
          ]
          if (state.var_qc_metrics) {
            newState += ["var_qc_metrics": state.var_qc_metrics]
          }
        return newState
        },
        toState: ["input": "output"]
      )
      | publish.run(
        fromState: { id, state -> [
            "input": state.input,
            "output": state.workflow_output,
            "compression": "gzip"
          ]
        },
        toState: ["output": "output"]
      )
      | setState(["output"]) 

  emit:
  output_ch
}
