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
          "layer": state.input_layer,
          "modality": state.modality,
          "var_name_filter": "scgpt_filter_with_hvg",
          "n_top_features": state.n_hvg,
          "flavor": state.hvg_flavor
        ]
      },
      toState: ["input": "output"]
    )
    | cross_check_genes.run(
      fromState: { id, state -> [
      // Check whether the genes are part of the provided vocabulary.
          "input": state.input,
          "modality": state.modality,
          "vocab_file": state.model_vocab,
          "input_var_gene_names": state.var_gene_names,
          "output": state.output,
          "pad_token": state.pad_token,
          "var_input": "scgpt_filter_with_hvg",
          "output_var_filter": "scgpt_cross_checked_genes"
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
          "output_obsm_binned_counts": "binned_counts",
          "var_input": "scgpt_cross_checked_genes",
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
          "input_obsm_binned_counts": "binned_counts",
          "var_input": "scgpt_cross_checked_genes",
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
          "obsm_embeddings": state.obsm_integrated,
          "finetuned_checkpoints_key": state.finetuned_checkpoints_key,
          "output": state.output
        ]
      },
      toState: ["input": "output"]
    )
    
    | neighbors_leiden_umap.run(
      fromState: [
        "input": "input",
        "obsm_input": "obsm_integrated",
        "modality": "modality",
        "uns_neighbors": "uns_neighbors",
        "obsp_neighbor_distances": "obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "obsp_neighbor_connectivities",
        "output": "workflow_output",
        "leiden_resolution": "leiden_resolution",
        "obsm_umap": "obsm_integrated",
      ],
      toState: ["output": "output"],
      args: [
        "uns_neighbors": "scGPT_integration_neighbors",
        "obsp_neighbor_distances": "scGPT_integration_distances",
        "obsp_neighbor_connectivities": "scGPT_integration_connectivities",
        "obs_cluster": "scGPT_integration_leiden",
        "obsm_umap": "X_scGPT_umap",
      ]
    )
    | setState(["output"])
  
  emit:
    output_ch
}