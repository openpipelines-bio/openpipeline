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
    | highly_variable_features_scanpy.run(
      fromState: {id, state ->
      // Annotates the mudata object with highly variable genes.
        [
          "input": state.input,
          "layer": state.layer,
          "modality": state.modality,
          "var_name_filter": "filter_with_hvg",
          "n_top_features": state.n_hvg,
          "flavor": "seurat_v3"
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      fromState: {id, state ->
        // do_filter does not need a layer argument because it filters all layers
        // from a modality.
        // filters the mudata object based on the HVG
        [
          "input": state.input,
          "modality": state.modality,
          "var_filter": "filter_with_hvg"
        ]
      },
      toState: ["input": "output"]
    )
    | cross_check_genes.run(
      fromState: { id, state ->
      // Check whether the genes are part of the provided vocabulary. Subsets for genes present in vocab only.
        [
          "input": state.input,
          "modality": state.modality,
          "vocab_file": state.model_vocab,
          "var_gene_names": state.var_gene_names,
          "output": state.output,
          "pad_token": state.pad_token
        ]
      },
      toState: ["input": "output"]
    )
    | binning.run(
      // Bins the data into a fixed number of bins.
        fromState: {id, state -> [
            "input": state.input,
            "modality": state.modality,
            "input_layer": state.input_layer,
            "n_input_bins": state.n_input_bins,
            "binned_layer": "binned",
            "output": state.output
          ]
        },
        toState: ["input": "output"]
    )
    | pad_tokenize.run(
      // Padding and tokenization of gene count values.
       fromState: {id, state -> [
            "input": state.input,
            "modality": state.modality,
            "model_vocab": state.model_vocab,
            "input_layer": "binned",
            "var_gene_names": state.var_gene_names,
            "pad_token": state.pad_token,
            "pad_value": state.pad_value,
            "max_seq_len": state.max_seq_len,
            "obsm_gene_tokens": "gene_id_tokens",
            "obsm_tokenized_values": "values_tokenized",
            "obsm_padding_mask": "padding_mask",
            "output": state.output
          ]
        },
        toState: ["input": "output"]
    )
    | embedding.run(
      // Generation of cell embedings from the tokenized gene counts values.
      fromState: {id, state -> [
          "input": state.input,
          "modality": state.modality,
          "model": state.model,
          "model_vocab": state.model_vocab,
          "model_config": state.model_config,
          "obsm_gene_tokens": "gene_id_tokens",
          "obsm_tokenized_values": "values_tokenized",
          "obsm_padding_mask": "padding_mask",
          "var_gene_names": state.var_gene_names,
          "obs_batch_label": state.obs_batch_label,
          "pad_token": state.pad_token,
          "pad_value": state.pad_value,
          "dsbn": state.dsbn,
          "batch_size": state.batch_size,
          "obsm_embeddings": "X_scGPT",
          "output": state.output
        ]
      },
      toState: ["input": "output"]
    )

    | find_neighbors.run(
      fromState: {id, state -> [
          "input": state.input,
          "uns_output": "scGPT_integration_neighbors",
          "obsp_distances": "scGPT_integration_distances",
          "obsp_connectivities": "scGPT_integration_connectivities",
          "obsm_input": "X_scGPT",
          "modality": state.modality
        ]
      },
      toState: ["input": "output"]
    )

    | leiden.run(
      runIf: {id, state -> state.var_name_mitochondrial_genes}
      fromState: {id, state -> [
        "input": state.input,
        "obsp_connectivities": "scGPT_integration_connectivities",
        "obsm_name": "scGPT_integration_leiden",
        "resolution": state.leiden_resolution,
        "modality": state.modality,
        ]
      },
      toState: ["input": "output"]
    )
    
    | move_obsm_to_obs.run(
      runIf: {id, state -> state.var_name_mitochondrial_genes}
      fromState: {id, state -> [
          "input": state.input,
          "obsm_key": "scGPT_integration_leiden",
          "modality": state.modality,
        ]
      },
      toState: ["input": "output"]
    )

    | umap.run(
      fromState: {id, state -> [
          "input": state.input,
          "uns_neighbors": "scGPT_integration_neighbors",
          "obsm_output": "X_scGPT_umap",
          "modality": state.modality,
          "output_compression": state.output_compression,
          "output": state.workflow_output
        ]
      },
      toState: { id, output, state ->
        [ output: output.output ]
      },
      auto: [ publish: true ]
    )
  
  emit:
    output_ch
}
