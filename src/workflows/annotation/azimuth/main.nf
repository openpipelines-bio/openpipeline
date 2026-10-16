workflow run_wf {
  take:
    input_ch

  main:
    prep_ch = input_ch
      | map {id, state ->
        assert state.input || state.input_seurat :
          "Either --input (h5mu) or --input_seurat (rds) must be provided."
        assert !(state.input && state.input_seurat) :
          "Provide only one of --input or --input_seurat."
        def n_levels = state.annotation_levels.size()
        assert !state.output_obs_predictions || state.output_obs_predictions.size() == n_levels :
          "--output_obs_predictions should have the same length as --annotation_levels (${n_levels})."
        assert !state.output_obs_probability || state.output_obs_probability.size() == n_levels :
          "--output_obs_probability should have the same length as --annotation_levels (${n_levels})."
        [id, state]
      }
      | map {id, state ->
        def from_h5mu = state.input != null
        // Viash fills in a '$id.$key.output_seurat.rds' template for optional outputs
        // that were not provided. Only keep --output_seurat if it was requested.
        def output_seurat = state.output_seurat
        if (output_seurat && output_seurat.toString().contains('$key')) {
          output_seurat = null
        }

        def seurat_layer = from_h5mu ? (state.input_layer ?: "X") : (state.input_layer ?: "counts")
        def new_state = state + [
          // Set aside the output for this workflow to avoid conflicts
          "workflow_output": state.output,
          "output_seurat": output_seurat,
          "seurat_assay": from_h5mu ? "RNA" : state.input_seurat_assay,
          "seurat_layer": seurat_layer,
          "seurat_x_layer": from_h5mu ? "X" : seurat_layer
        ]
        [id, new_state]
      }

    // H5MU input: split the modalities so that only the selected modality
    // is converted and annotated, the other modalities are merged back later.
    modalities_ch = prep_ch
      | filter {id, state -> state.input}
      | split_modalities.run(
        fromState: {id, state -> ["input": state.input, "id": id]},
        toState: ["split_output": "output", "split_output_types": "output_types"]
      )
      | flatMap {id, state ->
        def csv = state.split_output_types.splitCsv(strip: true, sep: ",").findAll{row -> !row[0].startsWith("#")}
        def header = csv.head()
        def types = csv.tail().collect { row ->
          [header, row].transpose().collectEntries()
        }
        def found_modalities = types.collect{dat -> dat.name}
        assert state.modality in found_modalities :
          "Modality '${state.modality}' was not found in --input. Found modalities: ${found_modalities}."
        types.collect{ dat ->
          // The id is kept the same for all modalities, the events are grouped back together by id before merging
          def new_state = state + ["input": state.split_output.resolve(dat.filename), "current_modality": dat.name]
          [id, new_state]
        }
      }
      // Remove arguments from split modalities from state
      | map {id, state ->
        def keysToRemove = ["split_output", "split_output_types"]
        [id, state.findAll{k, v -> k !in keysToRemove}]
      }

    h5mu_query_ch = modalities_ch
      | filter {id, state -> state.current_modality == state.modality}
      | from_h5mu_or_h5ad_to_seurat.run(
        fromState: {id, state -> [
          "input": state.input,
          "modality": state.modality,
          "assay": state.seurat_assay,
          "x_mapping": state.seurat_x_layer
        ]},
        toState: ["input_seurat": "output"]
      )

    other_mod_ch = modalities_ch
      | filter {id, state -> state.current_modality != state.modality}

    // Seurat input: nothing to split, the object is annotated directly
    seurat_query_ch = prep_ch
      | filter {id, state -> state.input_seurat}
      | map {id, state -> [id, state + ["current_modality": state.modality]]}

    annotated_ch = h5mu_query_ch.mix(seurat_query_ch)
      | azimuth_seurat.run(
        fromState: {id, state ->
          def args = [
            "input": state.input_seurat,
            "assay": state.seurat_assay,
            "input_layer": state.seurat_layer,
            "input_var_gene_names": state.input_var_gene_names,
            "reference": state.reference,
            "annotation_levels": state.annotation_levels,
            "do_adt": state.do_adt,
            "k_weight": state.k_weight,
            "n_trees": state.n_trees,
            "mapping_score_k": state.mapping_score_k,
            "output_obs_predictions": state.output_obs_predictions,
            "output_obs_probability": state.output_obs_probability,
            "output_obs_mapping_score": state.output_obs_mapping_score,
            "output_obsm_umap": state.output_obsm_umap,
            "output_obsm_integrated": state.output_obsm_integrated
          ]
          if (state.output_seurat) {
            args.output = state.output_seurat
          }
          args
        },
        toState: ["azimuth_output": "output"]
      )

    annotated_mod_ch = annotated_ch
      | from_seurat_to_h5mu.run(
        fromState: {id, state -> [
          "input": state.azimuth_output,
          "assay": state.seurat_assay,
          "modality": state.modality,
          "x_mapping": state.seurat_x_layer
        ]},
        toState: ["input": "output"]
      )
      // The layer mapped to .X is also copied to .layers, remove the duplicate
      | delete_layer.run(
        fromState: {id, state -> [
          "input": state.input,
          "modality": state.modality,
          "layer": [state.seurat_x_layer],
          "missing_ok": true
        ]},
        toState: ["input": "output"]
      )

    // The imputed protein expression (--do_adt) is stored as a separate modality.
    imputed_adt_ch = annotated_ch
      | filter {id, state -> state.do_adt}
      | from_seurat_to_h5mu.run(
        key: "from_seurat_to_h5mu_imputed_adt",
        fromState: {id, state -> [
          "input": state.azimuth_output,
          "assay": "impADT",
          "modality": state.output_modality_imputed_adt,
          "x_mapping": "data",
          // Cell metadata and embeddings are already stored in the annotated modality
          "layers_mapping": "False",
          "obs_mapping": "False",
          "obsm_mapping": "False",
          "obsp_mapping": "False",
          "varm_mapping": "False",
          "varp_mapping": "False",
          "uns_mapping": "False"
        ]},
        toState: {id, output, state ->
          state + ["input": output.output, "current_modality": state.output_modality_imputed_adt]
        }
      )

    output_ch = annotated_mod_ch.mix(imputed_adt_ch, other_mod_ch)
      | groupTuple(by: 0, sort: "hash")
      | map {id, states ->
        def modalities = states.collect{s -> s.current_modality}
        assert modalities.size() == modalities.unique(false).size() :
          "Found duplicate modality names: ${modalities}. Make sure that --output_modality_imputed_adt does not already exist in --input."
        // All states are copies of the same state, except for the input and current_modality keys
        def annotated_state = states.find{s -> s.current_modality == s.modality}
        [id, annotated_state + ["input": states.collect{s -> s.input}]]
      }
      | merge.run(
        fromState: [
          "input": "input",
          "output": "workflow_output",
          "output_compression": "output_compression"
        ],
        toState: ["output": "output"]
      )
      | map {id, state ->
        def new_state = state + ["output_seurat": state.output_seurat ? state.azimuth_output : null]
        [id, new_state]
      }
      | setState(["output", "output_seurat"])

  emit:
    output_ch
}
