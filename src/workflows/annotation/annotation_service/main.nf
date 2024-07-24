workflow process_reference {
  take:
    input_ch

  main:
    reference_ch = input_ch

    // Create reference specific output for this channel
    | map {id, state ->
      def new_state = state + ["reference_processed": state.output]
      [id, new_state]
      }

    // convert the reference h5ad file to h5mu
    | from_h5ad_to_h5mu.run(
        fromState: { id, state ->
        [
          "input": state.reference,
          "modality": "rna",
        ]
      },
      toState: [
        "input": "output",
      ]
    )
    // Split reference based on batches
    | split_samples.run(
        fromState: { id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "obs_feature": state.obs_reference_batch
        ]
      },
      toState: [ "output": "output", "output_files": "output_files" ]
      ) 
    // Turn each batch input h5mu into channel event
    | flatMap {id, state ->
        def outputDir = state.output
        def files = readCsv(state.output_files.toUriString())
        files.collect{ dat ->
          def new_id = dat.name
          def new_data = outputDir.resolve(dat.filename)
          [ new_id, state + ["reference_input": new_data]]
        }
      }
    // run process_samples workflow on reference and publish the processed reference
    | process_samples_workflow.run(
      fromState: {id, state ->
        def newState = [
          "input": state.reference_input, 
          "id": id,
          "rna_layer": state.reference_rna_layer,
          "add_id_to_obs": "true",
          "add_id_obs_output": "sample_id",
          "add_id_make_observation_keys_unique": "false",
          "rna_min_counts": state.rna_min_counts,
          "rna_max_counts": state.rna_max_counts,
          "rna_min_genes_per_cell": state.rna_min_genes_per_cell,
          "rna_max_genes_per_cell": state.rna_max_genes_per_cell,
          "rna_min_cells_per_gene": state.rna_min_cells_per_gene,
          "rna_min_fraction_mito": state.rna_min_fraction_mito,
          "rna_max_fraction_mito": state.rna_max_fraction_mito,
          "highly_variable_features_var_output": state.highly_variable_features_var_output,
          "highly_variable_features_obs_batch_key": state.highly_variable_features_obs_batch_key,
          "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
          "var_gene_names": state.var_gene_names,
          "mitochondrial_gene_regex": state.mitochondrial_gene_regex,
          "var_qc_metrics": state.var_qc_metrics,
          "top_n_vars": state.top_n_vars,
          "pca_overwrite": "true"
          ]
      },
      toState: {id, output, state -> ["reference_processed": output.output]},
      )

  emit:
    reference_ch
}

workflow process_query {
  take:
    input_ch

  main:
    query_ch = input_ch
    // Create query specific output for this channel
    | map {id, state ->
      def new_state = state + ["query_processed": state.output]
      [id, new_state]
      }
    // consider individual input files as events for process_samples pipeline
    | flatMap {id, state ->
        def outputDir = state.query_output
        def query_files = state.input
        // Workflow can take multiple input files. Split into seperate events.
        query_files.collect{ dat ->
          def filename = dat.getName()
          // make id's unique based on filename
          def new_id = filename.substring(0, filename.lastIndexOf('.h5mu'))
          def new_data = dat
          [ new_id, state + ["query_input": new_data]]
        }
      }
    | process_samples_workflow.run(
      fromState: {id, state ->
        def newState = [
          "input": state.query_input, 
          "id": id,
          "rna_layer": state.query_rna_layer,
          "add_id_to_obs": "true",
          "add_id_obs_output": "sample_id",
          "add_id_make_observation_keys_unique": "false",
          "rna_min_counts": state.rna_min_counts,
          "rna_max_counts": state.rna_max_counts,
          "rna_min_genes_per_cell": state.rna_min_genes_per_cell,
          "rna_max_genes_per_cell": state.rna_max_genes_per_cell,
          "rna_min_cells_per_gene": state.rna_min_cells_per_gene,
          "rna_min_fraction_mito": state.rna_min_fraction_mito,
          "rna_max_fraction_mito": state.rna_max_fraction_mito,
          "highly_variable_features_var_output": state.highly_variable_features_var_output,
          "highly_variable_features_obs_batch_key": state.highly_variable_features_obs_batch_key,
          "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
          "var_gene_names": state.var_gene_names,
          "mitochondrial_gene_regex": state.mitochondrial_gene_regex,
          "var_qc_metrics": state.var_qc_metrics,
          "top_n_vars": state.top_n_vars,
          "pca_overwrite": "true"
          ]
      },
      toState: {id, output, state -> ["query_processed": output.output]}, 
      )
    
  emit:
    query_ch
}

workflow run_wf {
  take:
    input_ch


  main:
    reference_ch = input_ch
      | process_reference
      | view {"After processing reference: $it"}

    query_ch = input_ch
      | process_query
      | view {"After processing query: $it"}

    // add id as _meta join id to be able to merge with source channel and end of workflow
    input_id_ch = input_ch
      | map{ id, state -> 
        def new_state = state + ["_meta": ["join_id": id]]
        [id, new_state]
        }

    // Mix the input channel with the processed query
    processed_ch = input_id_ch.mix(query_ch).mix(reference_ch)
  
    output_ch = processed_ch
      | view {"After processing: $it"}
      // Create workflow specific output for this channel
      // Make sure that process_query and process_reference have same output id
      | map { id, state -> ["processed", state]}
      // Combine output states of process_query and process_reference
      | groupTuple(by: 0)
      | map { id, state -> 
        def newState = state.collectEntries{it}
          [id, newState]
        }
      | view {"After mapping: $it"}
      | map {id, state ->
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
        }
      | view {"After mixing: $it"}
      | harmony_knn_workflow.run(
        runIf: { id, state -> state.annotation_methods.contains("harmony_knn") },
        fromState: { id, state ->
          def output_obs_predictions = state.obs_reference_targets.collect{it + "_pred_knn_harmony"}
          def output_obs_probabilities = state.obs_reference_targets.collect{it + "_proba_knn_harmony"}
          [ 
            "id": id,
            "input_query_dataset": state.query_processed,
            "input_reference_dataset": state.reference_processed,
            "modality": "rna",
            "embedding": "X_pca",
            "obs_reference_targets": state.obs_reference_targets,
            "output_obs_predictions": output_obs_predictions,
            "output_obs_probability": output_obs_probabilities,
            "output_compression": state.output_compression,
            "theta": state.theta,
            "obsm_integrated": "X_integrated_harmony",
            "obs_covariates": ["sample_id"],
            "weights": state.weights,
            "n_neighbors": state.n_neighbors,
            "output": state.output
          ]
        },
        toState: [ "query_processed": "output" ]
        )
      | scgpt_annotation_workflow.run(
        runIf: { id, state -> state.annotation_methods.contains("scgpt_annotation") },
        fromState: { id, state ->
          [ 
            "id": id,
            "input": state.query_processed,
            "modality": "rna",
            "input_layer": state.query_rna_layer,
            "var_gene_names": state.var_query_gene_names,
            "obs_batch_label": state.obs_query_batch,
            "model": state.model,
            "model_config": state.model_config,
            "model_vocab": state.model_vocab,
            "finetuned_checkpoints_key": state.finetuned_checkpoints_key,
            "label_mapper_key": state.label_mapper_key,
            "output_obs_predictions": "scgpt_pred",
            "output_obs_probability": "scgpt_proba",
            "pad_token": state.pad_token,
            "pad_value": state.pad_value,
            "n_hvg": state.n_hvg,
            "dsbn": state.dsbn,
            "batch_size": state.batch_size,
            "n_input_bins": state.n_input_bins,
            "seed": state.seed
          ]
        },
        toState: [ "query_processed": "output" ]
        )
        | scgpt_integration_knn_workflow.run(
        runIf: { id, state -> state.annotation_methods.contains("scgpt_annotation") },
        fromState: { id, state ->
          [ 
            "id": id,
            "input": state.query_processed,
            "modality": "rna",
            "input_layer": state.query_rna_layer,
            "var_gene_names": state.var_query_gene_names,
            "obs_batch_label": state.obs_query_batch,
            "model": state.model,
            "model_config": state.model_config,
            "model_vocab": state.model_vocab,
            "finetuned_checkpoints_key": state.finetuned_checkpoints_key,
            "label_mapper_key": state.label_mapper_key,
            "output_obs_predictions": "scgpt_pred",
            "output_obs_probability": "scgpt_proba",
            "pad_token": state.pad_token,
            "pad_value": state.pad_value,
            "n_hvg": state.n_hvg,
            "dsbn": state.dsbn,
            "batch_size": state.batch_size,
            "n_input_bins": state.n_input_bins,
            "seed": state.seed
          ]
        },
        toState: [ "query_processed": "output" ]
        )
      | setState(["output": "query_processed", "_meta": "_meta"])

  emit:
    output_ch
}
