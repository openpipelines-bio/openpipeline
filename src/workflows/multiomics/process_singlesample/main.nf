workflow run_wf {
  take:
    input_ch

  main:
    modalities_ch = input_ch
      // Make sure there is not conflict between the output from this workflow
      // And the output from any of the components
      | map {id, state ->
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
      }
      // If requested, add the id of the events (samples) a column in .obs. 
      // Also allows to make .obs_names (the .obs index) unique, by prefixing the values with an unique id per .h5mu file.
      // The latter is usefull to avoid duplicate observations during concatenation.
      | add_id.run(
        filter: {id, state -> state.add_id_to_obs },
        fromState: {id, state -> 
          def newState = [
            "input": state.input,
            "input_id": id,
            "make_observation_keys_unique": state.add_id_make_observation_keys_unique,
            "obs_output": state.add_id_obs_output,
            "add_id_to_obs": state.add_id_to_obs
          ]
          newState
        },
        toState: {id, output, state ->
          def keysToRemove = ["add_id_to_obs", "add_id_obs_output", "add_id_make_observation_keys_unique"]
          def newState = state.findAll{it.key !in keysToRemove}
          newState + ["input": output.output]
        }
      )
      | split_modalities_workflow.run(
        fromState: {id, state ->
          def newState = ["input": state.input, "id": id]
        },
        toState: ["output": "output", "output_types": "output_types"]
      )
      | flatMap {id, state ->
        def outputDir = state.output
        def csv = state.output_types.splitCsv(strip: true, sep: ",").findAll{!it[0].startsWith("#")}
        def header = csv.head()
        def types = csv.tail().collect { row ->
            [header, row].transpose().collectEntries()
        }
        
        types.collect{ dat ->
          // def new_id = id + "_" + dat.name
          def new_id = id // it's okay because the channel will get split up anyways
          def new_data = outputDir.resolve(dat.filename)
          [ new_id, state + ["input": new_data, modality: dat.name]]
        }
      }
      // Remove arguments from split modalities from state
      | map {id, state -> 
        def keysToRemove = ["output_types"]
        def newState = state.findAll{it.key !in keysToRemove}
        [id, newState]
      }
      | view {"After splitting modalities: $it"}


    //
    // Singlesample processing
    //
    def singlesample_arguments = [
      "rna": [
        "min_counts": "rna_min_counts",
        "max_counts": "rna_max_counts",
        "min_genes_per_cell": "rna_min_genes_per_cell",
        "max_genes_per_cell": "rna_max_genes_per_cell",
        "min_cells_per_gene": "rna_min_cells_per_gene",
        "min_fraction_mito": "rna_min_fraction_mito",
        "max_fraction_mito": "rna_max_fraction_mito",
        "var_name_mitochondrial_genes": "var_name_mitochondrial_genes",
        "obs_name_mitochondrial_fraction": "obs_name_mitochondrial_fraction",
        "var_gene_names": "var_gene_names",
        "mitochondrial_gene_regex": "mitochondrial_gene_regex",
        "min_fraction_ribo": "rna_min_fraction_ribo",
        "max_fraction_ribo": "rna_max_fraction_ribo",
        "var_name_ribosomal_genes": "var_name_ribosomal_genes",
        "obs_name_ribosomal_fraction": "obs_name_ribosomal_fraction",
        "ribosomal_gene_regex": "ribosomal_gene_regex",
        "layer": "rna_layer",
        "skip_scrublet_doublet_detection": "skip_scrublet_doublet_detection"
      ],
      "prot": [
        "min_counts": "prot_min_counts",
        "max_counts": "prot_max_counts",
        "min_proteins_per_cell": "prot_min_proteins_per_cell",
        "max_proteins_per_cell": "prot_max_proteins_per_cell",
        "min_cells_per_protein": "prot_min_cells_per_protein",
        "layer": "prot_layer",
      ],
      "gdo": [
        "min_counts": "gdo_min_counts",
        "max_counts": "gdo_max_counts",
        "min_guides_per_cell": "gdo_min_guides_per_cell",
        "max_guides_per_cell": "gdo_max_guides_per_cell",
        "min_cells_per_guide": "gdo_min_cells_per_guide",
        "layer": "gdo_layer",
      ], 
    ].asImmutable()

    multisample_ch_known = modalities_ch 
      // run the singlesample processing
      | runEach(
        components: [rna_singlesample, prot_singlesample, gdo_singlesample],
        filter: { id, state, component ->
          state.modality + "_singlesample" == component.config.name
        },
        fromState: { id, state, component ->
          def newState = singlesample_arguments.get(state.modality).collectEntries{key_, value_ -> 
            [key_, state[value_]]
          }
          return newState + ["id": id, "input": state.input]
        },
        toState: ["input": "output"],
      )

    multisample_ch_unknown = modalities_ch
      | filter{id, state -> state.modality !in singlesample_arguments.keySet()}

    output_ch = multisample_ch_unknown.mix(multisample_ch_known)
      // Group the modalities back together
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
        fromState: [
          "input": "input",
          "output": "workflow_output"
        ],
        toState: {id, output, state -> 
          ["output": output.output]
        }
      )

      | view {"After singlesample processing: $it"}
    // output_ch = modalities_ch

  emit:
    output_ch
}
