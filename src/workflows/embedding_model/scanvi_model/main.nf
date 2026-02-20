workflow run_wf {
  take:
    input_ch

  main:

    output_ch = input_ch
        // Set aside the output for this workflow to avoid conflicts
        | map {id, state -> 
        def new_state = state + ["workflow_output": state.output, "workflow_output_model": state.output_model]
        [id, new_state]
        }

        // Integrate & generate scvi model from the reference data
        | scvi.run(
          fromState: [
              "input": "input",
              "modality": "modality",
              "input_layer": "layer",
              "obs_batch": "obs_batch_label",
              "var_input": "var_hvg",
              "var_gene_names": "var_gene_names",
              "obs_size_factor": "obs_size_factor",
              "obs_categorical_covariate": "obs_categorical_covariate",
              "obs_continuous_covariate": "obs_continuous_covariate",
              "early_stopping": "early_stopping",
              "early_stopping_monitor": "early_stopping_monitor",
              "early_stopping_patience": "early_stopping_patience",
              "early_stopping_min_delta": "early_stopping_min_delta",
              "max_epochs": "max_epochs",
              "reduce_lr_on_plateau": "reduce_lr_on_plateau",
              "lr_factor": "lr_factor",
              "lr_patience": "lr_patience",
          ],
          args: [
              "obsm_output": "X_integrated_scvi"
          ],
          toState: [
              "reference": "output",
              "output_scvi_model": "output_model"
          ]
        )
        
        // Create scanvi model from the scvi reference model and integrate reference data
        | scanvi.run(
          fromState: [
              "input": "reference",
              "modality": "modality",
              "input_layer": "layer",
              "var_input": "var_hvg",
              "var_gene_names": "var_gene_names",
              "obs_labels": "obs_target",
              "unlabeled_category": "unlabeled_category",
              "scvi_model": "output_scvi_model",
              "early_stopping": "early_stopping",
              "early_stopping_monitor": "early_stopping_monitor",
              "early_stopping_patience": "early_stopping_patience",
              "early_stopping_min_delta": "early_stopping_min_delta",
              "max_epochs": "max_epochs",
              "reduce_lr_on_plateau": "reduce_lr_on_plateau",
              "lr_factor": "lr_factor",
              "lr_patience": "lr_patience",
          ],
          args: [
              "obsm_output": "X_integrated_scanvi",
              "obs_output_predictions": "scanvi_pred",
              "obs_output_probabilities": "scanvi_probabilities"
          ],
          toState: [
              "output": "output",
              "output_scanvi_model": "output_model"
          ]
        )

        | setState(["output", "output_scvi_model", "output_scanvi_model"])

  emit:
    output_ch
}