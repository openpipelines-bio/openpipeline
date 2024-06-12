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
    // Split reference based on batches
    | split_samples.run(
        fromState: { id, state ->
        [
          "input": state.reference,
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
          "output": "reference_processed.h5mu",
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
      toState: ["reference_processed": "output"],
      auto: [ publish: true ]
      )
    | map {id, state -> 
        def keysToRemove = ["output_files", "reference_input"]
        def newState = state.findAll{it.key !in keysToRemove}
        [id, newState]
      }
    | map {id, state -> 
        ["reference", state] 
      }

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
          "output": "query_processed.h5mu",
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
      toState: {id, output, state -> 
        ["query_processed": output.output]
      },
      auto: [ publish: true ]
      )
      | map { id, state -> 
        ["query", state] 
        }
    

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

    output_ch = reference_ch.mix(query_ch)
      | map { id, state -> ["processed", state]}
      | groupTuple(by: 0)
      | map { id, state -> [id, state.flatten()]}
      | niceView()


  emit:
    output_ch
}
