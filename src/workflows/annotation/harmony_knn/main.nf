workflow run_wf {
  take:
    input_ch

  main:

    // add id as _meta join id to be able to merge with source channel and end of workflow
    output_ch = input_ch
        // Set aside the output for this workflow to avoid conflicts
        | map {id, state -> 
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
        }
        // Add join_id to state
        | map{ id, state -> 
        def new_state = state + ["_meta": ["join_id": id]]
        [id, new_state]
        }
        | view {"After adding join_id: $it"}

        // Add 'query' id to .obs columns of query dataset
        | add_id.run(
            fromState: {id, state ->
                [
                "input": state.input,
                "input_id": "query",
                "obs_output": "dataset",
                ]
            },
            toState: ["input": "output"])
        // Add 'reference'id to .obs columns of reference dataset
        | add_id.run(
                fromState: {id, state ->
                    [
                    "input": state.reference,
                    "input_id": "reference",
                    "obs_output": "dataset",
                    ]
                },
                toState: ["reference": "output"])
        // Concatenate query and reference datasets
        | concatenate_h5mu.run(
            fromState: { id, state ->
            [
                "input": [state.input, state.reference],
                "input_id": ["query", "reference"],
                "other_axis_mode": "move"
            ]
            },
            toState: ["input": "output"]
            )
        | view {"After concatenation: $it"}
        // Run harmony integration with leiden clustering
        | harmony_leiden_workflow.run(
            fromState: { id, state ->
            [
                "id": id,
                "input": state.input,
                "modality": "rna",
                "uns_neighbors": "harmonypy_integration_neighbors",
                "obsp_neighbor_distances": "harmonypy_integration_distances",
                "obsp_neighbor_connectivities": "harmonypy_integration_connectivities",
                "embedding": state.obsm_embedding,
                "obsm_integrated": state.obsm_integrated,
                "theta": state.theta,
                "obs_covariates": state.obs_covariates,
                "obs_cluster": "harmony_integration_leiden",
                "leiden_resolution": state.leiden_resolution,
                "obsm_umap": "X_leiden_harmony_umap"
            ]
            },
            toState: ["input": "output"]
            )
        | view {"After integration: $it"}
        // Split integrated dataset back into a separate reference and query dataset
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
        | view {"After splitting query: $it"}
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