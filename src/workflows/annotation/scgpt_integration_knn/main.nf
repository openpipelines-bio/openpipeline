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
        // add id as _meta join id to be able to merge with source channel and end of workflow
        | map{ id, state -> 
        def new_state = state + ["_meta": ["join_id": id]]
        [id, new_state]
        }
        | view {"After adding join_id: $it"}
        // Add 'query' id to .obs columns of query dataset
        | add_id.run(
            fromState: [
                "input": "input",
            ],
            args:[
                "input_id": "query",
                "obs_output": "dataset",
            ],
            toState: ["input": "output"])
        // Add 'reference'id to .obs columns of reference dataset
        | add_id.run(
                fromState:[
                    "input": "reference",
                ],
                args:[
                    "input_id": "reference",
                    "obs_output": "dataset"
                ],
                toState: ["reference": "output"])
        // Make sure that query and reference dataset have batch information in the same .obs column
        // By copying the respective .obs columns to the obs column "batch_label"
        | copy_obs.run(
            fromState: [
                "input": "input",
                "modality": "modality",
                "input_obs_key": "input_obs_batch_label",
            ],
            args: [
                "output_obs_key": "batch_label"
            ],
            toState: [
                "input": "output"
            ]
        )
        | copy_obs.run(
            fromState: [
                "input": "reference",
                "modality": "modality",
                "input_obs_key": "reference_obs_batch_label",
            ],
            args: [
                "output_obs_key": "batch_label"
            ],
            toState: [
                "reference": "output"
            ]
        )
        // Make sure the .var columns containing the gene names has the same name prior to integration
        // By copying the respective .var columns to the var column "gene_symbols"
        | copy_var.run(
            fromState: [
                "input": "input",
                "modality": "modality",
                "input_var_key": "input_var_gene_names"
            ],
            args: [
                "output_var_key": "gene_symbols"
            ],
            toState: [
                "input": "output"
            ]
        )
        | copy_var.run(
            fromState: [
                "input": "reference",
                "modality": "modality",
                "input_var_key": "reference_var_gene_names"
            ],
            args: [
                "output_var_key": "gene_symbols"
            ],
            toState: [
                "reference": "output"
            ]
        )
        // Concatenate query and reference datasets prior to integration
        | concatenate_h5mu.run(
            fromState: { id, state -> [
                "input": [state.input, state.reference]
                ]
            },
            args: [
                "input_id": ["query", "reference"],
                "other_axis_mode": "move"
            ],
            toState: ["input": "output"]
            )
        | view {"After concatenation: $it"}
        // Run scgpt integration with leiden clustering
        | scgpt_leiden_workflow.run(
            fromState: { id, state ->
            [
                "id": id,
                "input": state.input,
                "modality": state.modality,
                "input_layer": state.input_layer,
                "model": state.model,
                "model_vocab": state.model_vocab,
                "model_config": state.model_config,
                "finetuned_checkpoints_key": state.finetuned_checkpoints_key,
                "obsm_integrated": state.output_obsm_integrated,
                "pad_token": state.pad_token,
                "pad_value": state.pad_value,
                "n_hvg": state.n_hvg,
                "max_seq_len": state.max_seq_len,
                "dsbn": state.dsbn,
                "batch_size": state.batch_size,
                "n_input_bins": state.n_input_bins,
                "seed": state.seed,
                "leiden_resolution": state.leiden_resolution,
            ]},
            args: [
                "obs_batch_label": "batch_label",
                "var_gene_names": "gene_symbols"
            ],
            toState: ["input": "output"]
            )
        | view {"After integration: $it"}
        // Split integrated dataset back into a separate reference and query dataset
        | split_h5mu.run(
            fromState: [
                "input": "input",
                "modality": "modality"
            ],
            args: [
                "obs_feature": "dataset",
                "output_files": "sample_files.csv",
                "drop_obs_nan": "true",
                "output": "ref_query"
            ],
            toState: [ 
                "output": "output", 
                "output_files": "output_files" 
            ],
            )
        | view {"After sample splitting: $it"}
        // map the integrated query and reference datasets back to the state
        | map {id, state ->
            def outputDir = state.output
            def files = readCsv(state.output_files.toUriString())
            def query_file = files.findAll{ dat -> dat.name == 'query' }
            assert query_file.size() == 1, 'there should only be one query file'
            def reference_file = files.findAll{ dat -> dat.name == 'reference' }
            assert reference_file.size() == 1, 'there should only be one reference file'
            def integrated_query = outputDir.resolve(query_file.filename)
            def integrated_reference = outputDir.resolve(reference_file.filename)
            def newKeys = ["integrated_query": integrated_query, "integrated_reference": integrated_reference]
            [id, state + newKeys]
            }
        | view {"After splitting query: $it"}
        // Perform KNN label transfer from reference to query
        | pynndescent_knn.run(
            fromState: [
                "input": "integrated_query",
                "modality": "modality",
                "input_obsm_features": "output_obsm_integrated",
                "reference": "integrated_reference",
                "reference_obsm_features": "output_obsm_integrated",
                "reference_obs_targets": "reference_obs_targets",
                "output_obs_predictions": "output_obs_predictions",
                "output_obs_probability": "output_obs_probability",
                "output_compression": "output_compression",
                "weights": "weights",
                "n_neighbors": "n_neighbors",
                "output": "workflow_output"
            ],
            toState: {id, output, state -> ["output": output.output]},
            )
    
  emit:
    output_ch
}
