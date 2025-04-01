workflow run_wf {
  take:
  input_ch

  main:
    multisample_ch = input_ch
      // Make sure there is not conflict between the output from this workflow
      // And the output from any of the components
      | map {id, state ->
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
      }
      // The input for this workflow can either be a list of unimodal files
      // or a single multimodal file. To destingish between the two, the files will be split either way.
      // For multiple unimodal files, the result before or after splitting is identical.
      // In both cases, this workflow requires split files.
      
      // Split must be called on each item of the input list, so split it into multiple events with unique ids
      // Unique ids are required to run a component
      | flatMap {id, state ->
        def newEvents = state.input.withIndex().collect{input_file, index -> 
          def newState = state + ["input": input_file, "original_id": id]
          ["${id}_${index}", newState]
        }
        newEvents
      }
      | split_modalities_workflow.run(
        fromState: {id, state ->
          [
            "input": state.input,
            "id": id
          ]
        },
        toState: [
          "output": "output", 
          "output_types": "output_types"
        ]
      )
      // gather the output from split_modalities_workflow
      // by reading the output csv (the csv contains 1 line per output file)
      | flatMap {id, state ->
        def outputDir = state.output
        def types = readCsv(state.output_types.toUriString())
        
        types.collect{ dat ->
          def new_id = state.original_id + "_${dat.name}" // Make a unique ID by appending the modality name.
          def new_data = outputDir.resolve(dat.filename)
          [ new_id, state + ["input": new_data, modality: dat.name]]
        }
      }
      // Remove arguments from split modalities from state
      | map {id, state -> 
        def keysToRemove = ["output_types", "original_id"]
        def newState = state.findAll{it.key !in keysToRemove}
        [id, newState]
      }
    
    multisample_ch
      | toSortedList()
      | map{all_input ->
        def ids = all_input.collect({it[0]})
        assert ids.clone().unique().size() == ids.size(): "Found duplicate modalities in the input."
      }

    //
    // Multisample processing
    //
    def multisample_arguments = [
      "rna": [
        "highly_variable_features_var_output": "highly_variable_features_var_output",
        "highly_variable_features_obs_batch_key": "highly_variable_features_obs_batch_key",
        "var_qc_metrics": "var_qc_metrics",
        "top_n_vars": "top_n_vars",
        "layer": "rna_layer",
        "enable_scaling": "rna_enable_scaling",
        "scaling_output_layer": "rna_scaling_output_layer",
        "scaling_max_value": "rna_scaling_max_value",
        "scaling_zero_center": "rna_scaling_zero_center",
      ],
      "prot": [
        "layer": "prot_layer",
        "clr_axis": "clr_axis",
      ]
    ].asImmutable()

    multimodal_ch_known = multisample_ch
      | runEach(
        components: [rna_multisample, prot_multisample],
        filter: { id, state, component ->
          state.modality + "_multisample" == component.config.name
        },
        fromState: { id, state, component -> 
          def newState = multisample_arguments.get(state.modality).collectEntries{key_, value_ -> 
            [key_, state[value_]]
          }
          newState + ["id": id, "input": state.input]
        },
        toState: {id, output, state, component ->
          def newState = state + ["input": output.output]
          return newState
        }
      )

    multimodal_ch_unknown = multisample_ch
      | filter { id, state -> state.modality !in multisample_arguments.keySet() }
    
    multimodal_ch = multimodal_ch_unknown.mix(multimodal_ch_known)
      // Remove arguments for multisample processing from state.
      | map {id, state -> 
        def keysToRemove = multisample_arguments.inject([]){currentKeys, modality, stateMapping -> 
          def newKeys = currentKeys + stateMapping.values()
          return newKeys
        } 
        keysToRemove -= ["rna_enable_scaling", "rna_scaling_output_layer"]
        def newState = state.findAll{it.key !in keysToRemove }
        [id, newState]
      }
      | view {"After multisample processing: $it"}

      //
      // Merging: joining the observations from all modalities together. Everything in 1 file.
      //

      // Set the original IDs back into place to use them in groupTuple
      | map {id, state ->
        def newEvent = [state.id, state]
        newEvent
      }
      // Group the modalities back together per input sample
      | groupTuple(by: 0, sort: "hash")
      | view {"After toSortedList: $it"}
      | map { id, states ->
          def new_input = states.collect{it.input}
          def modalities = states.collect{it.modality}.unique()
          def other_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
          }.minus(["output", "input", "modality"])
          def new_state = other_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "Arguments should be the same across modalities. Please report this \
                                                 as a bug. Argument name: $argument_name, \
                                                 argument value: $argument_values"
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
          }
          [id, new_state + ["input": new_input, "modalities": modalities]]
      }
      | view {"Input merge channel: $it"}
      | merge.run(
        fromState: ["input": "input"],
        toState: ["input": "output"],
      )
      | view {"After merging processing: $it"}

      // Processing of multi-modal multisample MuData files.
      // Performs calculations on samples that have *not* been integrated,
      // and can be considered a "no-integration" workflow.
      output_ch = [dimensionality_reduction_rna, dimensionality_reduction_scaling_rna, dimensionality_reduction_prot].inject(multimodal_ch){ channel_in, component ->
        channel_out_integrated = channel_in
          | component.run(
            runIf: {id, state ->
              def reg = state.rna_enable_scaling ? ~/^dimensionality_reduction_(scaling_)?/ : ~/^dimensionality_reduction_/
              def modality_to_check = component.name - reg
              state.modalities.contains(modality_to_check)
            },
            fromState: { id, state -> 
              def stateMappings = [
                "dimensionality_reduction_rna": 
                  [
                    "id": id,
                    "input": state.input,
                    "layer": "log_normalized",
                    "modality": "rna",
                    "var_pca_feature_selection": state.highly_variable_features_var_output, // run PCA on highly variable genes only
                    "pca_overwrite": state.pca_overwrite,
                    "output": state.workflow_output,
                  ],
                "dimensionality_reduction_scaling_rna":
                  [
                    "id": id,
                    "input": state.input,
                    "layer": state.rna_scaling_output_layer,
                    "modality": "rna",
                    "var_pca_feature_selection": state.highly_variable_features_var_output, // run PCA on highly variable genes only
                    "pca_overwrite": state.pca_overwrite,
                    // extra scaling args
                    "obsm_pca": state.rna_scaling_pca_obsm_output,
                    "pca_loadings_varm_output": state.rna_scaling_pca_loadings_varm_output,
                    "pca_variance_uns_output": state.rna_scaling_pca_variance_uns_output,
                    "pca_overwrite": state.pca_overwrite,
                    "obsm_umap": state.rna_scaling_umap_obsm_output,
                    "uns_neighbors": "neighbors_scaled",
                    "obsp_neighbor_connectivities": "connectivities_scaled",
                    "obsp_neighbor_distances": "distances_scaled",
                    "output": state.workflow_output,
                  ],
                "dimensionality_reduction_prot":
                  [
                    "id": id,
                    "input": state.input,
                    "layer": "clr",
                    "modality": "prot",
                    "pca_overwrite": state.pca_overwrite,
                    "output": state.workflow_output,
                  ]
              ]
              return stateMappings[component.name]
            },
            toState: ["input": "output"]
          )
      }
      // At the end of the reduce statement,
      // the `toState` closure put the output back into
      // into the 'input' slot
      | map {id, state ->
        [id, ["output": state.input]]
      }

  emit:
  output_ch
}