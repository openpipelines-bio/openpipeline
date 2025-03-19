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
    // Annotates the mudata object with highly variable genes.
    | highly_variable_features_scanpy.run(
      fromState: [
          "input": "input",
          "layer": "input_layer",
          "modality": "modality",
          "n_top_features": "n_hvg",
          "flavor": "hvg_flavor"
        ],
      args: ["var_name_filter": "scgpt_filter_with_hvg"],
      toState: ["input": "output"]
    )
    // Check whether the genes are part of the provided vocabulary.
    | cross_check_genes.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "vocab_file": "model_vocab",
        "input_var_gene_names": "var_gene_names",
        "output": "output",
        "pad_token": "pad_token"
      ],
      args: [
        "var_input": "scgpt_filter_with_hvg",
        "output_var_filter": "scgpt_cross_checked_genes"
      ],
      toState: [
        "input": "output"
      ]
    )
    // Bins the data into a fixed number of bins.
    | binning.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "input_layer",
        "n_input_bins": "n_input_bins",
        "output": "output"
      ],
      args: [
        "output_obsm_binned_counts": "binned_counts",
        "var_input": "scgpt_cross_checked_genes"
      ],
      toState: [
        "input": "output"
      ]
    )
    // Padding and tokenization of gene count values.
    | pad_tokenize.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "model_vocab": "model_vocab",
        "var_gene_names": "var_gene_names",
        "pad_token": "pad_token",
        "pad_value": "pad_value",
        "max_seq_len": "max_seq_len",
        "output": "output"
      ],
      args: [
        "input_obsm_binned_counts": "binned_counts",
        "var_input": "scgpt_cross_checked_genes",
        "obsm_gene_tokens": "gene_id_tokens",
        "obsm_tokenized_values": "values_tokenized",
        "obsm_padding_mask": "padding_mask"
      ],
      toState: [
        "input": "output"
      ]
    )
    // Generation of cell embedings from the tokenized gene counts values.
    | embedding.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "model": "model",
        "model_vocab": "model_vocab",
        "model_config": "model_config",
        "var_gene_names": "var_gene_names",
        "obs_batch_label": "obs_batch_label",
        "pad_token": "pad_token",
        "pad_value": "pad_value",
        "dsbn": "dsbn",
        "batch_size": "batch_size",
        "obsm_embeddings": "obsm_integrated",
        "finetuned_checkpoints_key": "finetuned_checkpoints_key",
        "output": "output"
      ],
      args: [
        "obsm_gene_tokens": "gene_id_tokens",
        "obsm_tokenized_values": "values_tokenized",
        "obsm_padding_mask": "padding_mask"
      ],
      toState: [
        "input": "output"
      ]
    )
    // Calculation of neighbors, leiden clustering and UMAP.
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
      toState: [
        "output": "output"
      ],
      args: [
        "uns_neighbors": "scGPT_integration_neighbors",
        "obsp_neighbor_distances": "scGPT_integration_distances",
        "obsp_neighbor_connectivities": "scGPT_integration_connectivities",
        "obs_cluster": "scGPT_integration_leiden",
        "obsm_umap": "X_scGPT_umap"
      ]
    )
    | setState(["output"])
  
  emit:
    output_ch
}