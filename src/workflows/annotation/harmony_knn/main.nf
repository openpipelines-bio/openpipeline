workflow query_wf {
take:
    input_ch

main:
    output_ch = input_ch
        // Add 'query' id to .obs columns
        | add_id.run(
            fromState: {id, state ->
                [
                "input": state.input_query_dataset,
                "input_id": "query",
                "obs_output": "dataset",
                ]
            },
            toState: ["input": "output"])
            // Keep only rna modality of reference, to avoid problems with concatenation
        | split_modalities_workflow.run(
            fromState: {id, state ->
            def newState = ["input": state.input, "id": id]
            },
            toState: ["query_dir": "output", "query_types": "output_types"]
            )
        | flatMap {id, state ->
            def outputDir = state.query_dir
            def types = readCsv(state.query_types.toUriString())
            types.findAll { dat -> dat.name == state.modality }.collect{ dat ->
            def new_data = outputDir.resolve(dat.filename)
            [ id, ["query": new_data]]
            }
            }

emit:
    output_ch
}

workflow reference_wf {
take:
    input_ch

main:
    output_ch = input_ch
        // Add 'reference'id to .obs columns
        | add_id.run(
            fromState: {id, state ->
                [
                "input": state.input_reference_dataset,
                "input_id": "reference",
                "obs_output": "dataset",
                ]
            },
            toState: ["input": "output"])
            // Keep only rna modality of reference, to avoid problems with concatenation
        | split_modalities_workflow.run(
            fromState: {id, state ->
            def newState = ["input": state.input, "id": id]
            },
            toState: ["reference_dir": "output", "reference_types": "output_types"]
            )
        | flatMap {id, state ->
            def outputDir = state.reference_dir
            def types = readCsv(state.reference_types.toUriString())
            types.findAll { dat -> dat.name == state.modality }.collect{ dat ->
            def new_data = outputDir.resolve(dat.filename)
            [ id, ["reference": new_data]]
            }
            }

emit:
    output_ch
}

workflow run_wf {
  take:
    input_ch

  main:
    query_ch = input_ch
        | query_wf
        | view {"After processing query: $it"}

    reference_ch = input_ch
        | reference_wf
        | view {"After processing reference: $it"}
    
    // add id as _meta join id to be able to merge with source channel and end of workflow
    input_id_ch = input_ch
      | map{ id, state -> 
        def new_state = state + ["_meta": ["join_id": id]]
        [id, new_state]
      }
      | view {"After adding join_id: $it"}

    // Join the input channel with the processed query to synchronize outputs of query_ch and reference_ch
    processed_ch = input_id_ch.join(query_ch).join(reference_ch)
        | flatMap { id, input_state, query_state, reference_state ->
            // def state = tuple[1]
            // println "STATE: $state"
            def new_state = input_state + query_state + reference_state
            println "NEW_STATE: $new_state"
            [[id, new_state]]
        }
        | view {"After joining: $it"}
        

    output_ch = processed_ch
        // Concatenate query and reference datasets
        | concatenate_h5mu.run(
            fromState: { id, state ->
            [
                "input": [state.query, state.reference],
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
