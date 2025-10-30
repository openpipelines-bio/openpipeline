workflow run_wf {
  take:
    input_ch

  main:
    
    query_ch = input_ch
      // Log normalize query dataset to target sum of 10000
      | normalize_total.run(
        fromState: { id, state -> [
          "input": state.input,
          "modality": state.modality,
          "input_layer": state.input_layer,
        ]},
        args: [
          "output_layer": "normalized_10k",
          "target_sum": "10000",
        ],
        toState: [
          "input": "output",
        ]
      )
      | log1p.run( 
        fromState: { id, state -> [
          "input": state.input,
          "modality": state.modality
        ]},
        args: [
          "input_layer": "normalized_10k",
          "output_layer": "log_normalized_10k",
        ],
        toState: [
          "input": "output"
        ]
      )
      | delete_layer.run(
        fromState: { id, state -> [
          "input": state.input,
          "modality": state.modality
        ]},
        args: [
          "layer": "normalized_10k"
        ],
        toState: [
          "input": "output"
        ]
      )
      | view {"After query normalization: $it"}

    ref_ch = input_ch
      // Log normalize reference dataset to target sum of 10000
      | normalize_total.run(
        key: "normalize_total_reference",
        runIf: { id, state ->
          state.reference
        },
        fromState: { id, state -> [
          "input": state.reference,
          "modality": state.modality,
          "input_layer": state.reference_layer,
        ]},
        args: [
          "output_layer": "normalized_10k",
          "target_sum": "10000",
        ],
        toState: [
          "reference": "output",
        ]
      )
      | log1p.run( 
        key: "log1p_reference",
        runIf: { id, state ->
          state.reference
        },
        fromState: { id, state -> [
          "input": state.reference,
          "modality": state.modality
        ]},
        args: [
          "input_layer": "normalized_10k",
          "output_layer": "log_normalized_10k",
        ],
        toState: [
          "reference": "output"
        ]
      )
      | view {"After reference normalization: $it"}


    output_ch = query_ch.join(ref_ch, failOnMismatch: true, failOnDuplicate: true)
        | view {"After channel mixing: $it"}
        // Set aside the output for this workflow to avoid conflicts
        | map {id, query_state, ref_state -> 
          def newState = query_state + ["reference": ref_state.reference]
          [id, newState]
        }        
        // Run harmony integration with leiden clustering
        | celltypist_component.run(
          fromState: { id, state -> [
            "input": state.input,
            "modality": state.modality,
            "input_var_gene_names": state.input_var_gene_names,
            "input_reference_gene_overlap": state.input_reference_gene_overlap,
            "reference": state.reference,
            "reference_obs_target": state.reference_obs_target,
            "reference_var_gene_names": state.reference_var_gene_names,
            "reference_var_input": state.reference_var_input,
            "model": state.model,
            "feature_selection": state.feature_selection,
            "majority_voting": state.majority_voting,
            "C": state.C,
            "max_iter": state.max_iter,
            "use_SGD": state.use_SGD,
            "min_prop": state.min_prop,
            "output": state.output,
            "output_obs_predictions": state.output_obs_predictions,
            "output_obs_probability": state.output_obs_probability,
            "sanitize_ensembl_ids": state.sanitize_ensembl_ids
          ]},
          args: [
            "input_layer": "log_normalized_10k",
            "reference_layer": "log_normalized_10k"
          ],
          toState: [
            "output": "output"
          ]
        )
        | view {"After annotation: $it"}
        | setState(["output"])

  emit:
    output_ch
}
