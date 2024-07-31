workflow run_wf {

  take:
    input_ch

  main:
    output_ch = input_ch
    // Set aside the output for this workflow to avoid conflicts
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | highly_variable_features_scanpy.run(
      runIf: {id, state -> !state.model && state.n_hvg},
      fromState: {id, state ->
      // Annotates the mudata object with highly variable genes.
        [
          "input": state.input,
          "layer": state.input_layer,
          "modality": state.modality,
          "var_name_filter": "filter_with_hvg",
          "n_top_features": state.n_hvg,
          "flavor": "seurat_v3"
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      runIf: {id, state -> !state.model && state.n_hvg},
      fromState: {id, state ->
        // do_filter does not need a layer argument because it filters all layers
        // from a modality.
        // filters the mudata object based on the HVG
        [
          "input": state.input,
          "modality": state.modality,
          "var_filter": "filter_with_hvg"
        ]
      },
      toState: ["input": "output"]
    )
      | celltypist_annotation.run(
      // Run celltypist annotation
       fromState: {id, state -> [
          "input": state.input,
          "modality": state.modality,
          "input_layer": state.input_layer,
          "var_query_gene_names": state.var_query_gene_names,
          "reference": state.reference,
          "reference_layer": state.reference_layer,
          "reference_obs_targets": state.reference_obs_targets,
          "check_expression": state.check_expression,
          "var_reference_gene_names": state.var_reference_gene_names,
          "model": state.model,
          "feature_selection": state.feature_selection,
          "majority_voting": state.majority_voting,
          "output_obs_predictions": state.output_obs_predictions,
          "output_obs_probability": state.output_obs_probability,
          "output": state.workflow_output,
          "output_compression": state.output_compression
        ]
        },
        toState: {id, output, state -> ["output": output.output]},
        auto: [ publish: true]
    )

  emit:
    output_ch
}