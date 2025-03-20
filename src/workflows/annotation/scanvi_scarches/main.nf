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

        | scvi.run(
            fromState: [
                "input": "reference",
                "modality": "modality",
                "input_layer": "layer",
                "obs_batch": "reference_obs_batch_label",
                "var_input": "reference_var_hvg",
                "var_gene_names": "reference_var_gene_names",
                "obs_size_factor": "reference_obs_size_factor",
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
                "obsm_output": "X_scvi_integrated"
            ],
            toState: [
                "reference": "output",
                "output_scvi_model": "output_model"
            ]
        )
        | view {"After scvi: $it"}
        | scanvi.run(
            fromState: [
                "input": "reference",
                "modality": "modality",
                "input_layer": "layer",
                "var_input": "reference_var_hvg",
                "var_gene_names": "reference_var_gene_names",
                "obs_labels": "reference_obs_target",
                "unlabeled_category": "unlabeled_category",
                "scvi_model": "output_scvi_model",
                "obsm_output": "output_obsm_integrated",
                "obs_output_predictions": "output_obs_predictions",
                "obs_output_probabilities": "output_obs_probability",
                "early_stopping": "early_stopping",
                "early_stopping_monitor": "early_stopping_monitor",
                "early_stopping_patience": "early_stopping_patience",
                "early_stopping_min_delta": "early_stopping_min_delta",
                "max_epochs": "max_epochs",
                "reduce_lr_on_plateau": "reduce_lr_on_plateau",
                "lr_factor": "lr_factor",
                "lr_patience": "lr_patience",
            ],
            toState: [
                "reference": "output",
                "scanvi_model": "output_model"
            ]
        )

        | view {"After scanvi: $it"}
        | scarches.run(
            fromState: [
                "input": "input",
                "layer": "layer",
                "modality": "modality",
                "input_obs_batch": "input_obs_batch_label",
                "input_var_gene_names": "input_var_gene_names",
                "input_obs_size_factor": "input_obs_size_factor",
                "reference": "scanvi_model",
                "obsm_output": "output_obsm_integrated",
                "obs_output_predictions": "output_obs_predictions",
                "obs_output_probabilities": "output_obs_probability",
                "early_stopping": "early_stopping",
                "early_stopping_monitor": "early_stopping_monitor",
                "early_stopping_patience": "early_stopping_patience",
                "early_stopping_min_delta": "early_stopping_min_delta",
                "max_epochs": "max_epochs",
                "reduce_lr_on_plateau": "reduce_lr_on_plateau",
                "lr_factor": "lr_factor",
                "lr_patience": "lr_patience",
                "output": "workflow_output",
                "model_output": "workflow_output_model"
            ],
            toState: {id, output, state -> 
                [
                "output": output.output,
                "output_model": output.model_output
                ]
            }

        )
        | niceView()

  emit:
    output_ch
}