workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
        // Create workflow specific output for this channel
        | map {id, state ->
            def new_state = state + ["workflow_output": state.output]
            [id, new_state]
            }
        // Add 'reference' and 'query' id to .obs columns
        | add_id.run(
            fromState: {id, state ->
                [
                "input": state.input_query_dataset,
                "input_id": "query",
                "obs_output": "dataset",
                ]
            },
            toState: ["input_query_dataset": "output"])
        | add_id.run(
            fromState: {id, state ->
                [
                "input": state.input_reference_dataset,
                "input_id": "reference",
                "obs_output": "dataset",
                ]
            },
            toState: ["input_reference_dataset": "output"])
        // Concatenate query and reference datasets
        | concatenate_h5mu.run(
            fromState: { id, state ->
            [
                "input": [state.input_query_dataset, state.input_reference_dataset],
                "input_id": ["query", "reference"],
                "other_axis_mode": "move"
            ]
            },
            toState: ["input": "output"]
            )
        | view {"After concatenation: $it"}
        // Run harmony integration 
        | harmonypy.run(
            fromState: { id, state ->
            [
                "input": state.input,
                "modality": "rna",
                "obsm_input": state.embedding,
                "obsm_output": state.obsm_integrated,
                "theta": state.theta,
                "obs_covariates": state.obs_covariates
            ]
            },
            toState: ["input": "output"]
            )
        | view {"After integration: $it"}
        // Split integrated dataset back into reference and query dataset
        | split_samples.run(
            fromState: { id, state ->
            [
                "input": state.input,
                "modality": "rna",
                "obs_feature": "dataset",
                "output_files": "sample_files.csv",
                "drop_obs_nan": "true",
                "output": "ref_query"
            ]
            },
            toState: [ "output": "output", "output_files": "output_files" ],
            auto: [ publish: true ]
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
        | view {"After taking query: $it"}
        // Perform KNN label transfer from reference to query
        | pynndescent_knn.run(
            fromState: { id, state ->
            [
                "input": state.integrated_query,
                "modality": "rna",
                "input_obsm_features": state.obsm_integrated,
                "reference": state.integrated_reference,
                "reference_obsm_features": state.obsm_integrated,
                "reference_obs_targets": state.obs_reference_targets,
                "output_obs_predictions": state.output_obs_predictions,
                "output_obs_probability": state.output_obs_probability,
                "output_uns_parameters": state.output_uns_parameters,
                "output_compression": state.output_compression,
                "weights": state.weights,
                "n_neighbors": state.n_neighbors,
                "output": state.workflow_output
            ]
            },
            toState: {id, output, state -> ["output": output.output]},
            auto: [ publish: true ]
            )
    
  emit:
    output_ch
}
