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
          "input_var_gene_names": state.var_gene_names,
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
            "output_compression": state.output_compression,
            "binned_layer": state.binned_layer,
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
            "input_layer": state.binned_layer,
            "input_var_gene_names": state.var_gene_names,
            "n_input_bins": state.n_input_bins,
            "pad_token": state.pad_token,
            "pad_value": state.pad_value,
            "max_seq_len": state.max_seq_len,
            "output_compression": state.output_compression,
            "output_obsm_gene_tokens": state.obsm_gene_tokens,
            "output_obsm_tokenized_values": state.obsm_tokenized_values,
            "output_obsm_padding_mask": state.obsm_padding_mask,
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
          "input_obsm_gene_tokens": state.obsm_gene_tokens,
          "input_obsm_tokenized_values": state.obsm_tokenized_values,
          "input_obsm_padding_mask": state.obsm_padding_mask,
          "input_var_gene_names": state.var_gene_names,
          "input_obs_batch_id": state.obs_batch_id,
          "output": "workflow_output",
          "output_compression": state.output_compression,
          "embedding_layer_key": state.embedding_layer_key,
          "pad_token": state.pad_token,
          "pad_value": state.pad_value,
          "dropout": state.dropout,
          "DSBN": state.DSBN,
          "batch_size": state.batch_size,
        ]
      },
        auto: [ publish: true ]
    )
  
  emit:
    output_ch
}
