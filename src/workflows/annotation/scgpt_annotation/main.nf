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
    // Annotate the mudata object with highly variable genes.
    | highly_variable_features_scanpy.run(
      fromState: [
          "input": "input",
          "layer": "input_layer",
          "modality": "modality",
          "n_top_features": "n_hvg",
        ],
        args: [
          "var_name_filter": "filter_with_hvg",
          "flavor": "seurat_v3"
        ],
      toState: ["input": "output"]
    )
    // Filter the mudata object based on the HVG
    // do_filter does not need a layer argument because it filters all layers from a modality.
    | do_filter.run(
      fromState: [
          "input": "input",
          "modality": "modality"
        ],
        args:[
          "var_filter": "filter_with_hvg"
        ],
      toState: ["input": "output"]
    )
    // Check whether the genes are part of the provided vocabulary. 
    // Subsets for genes present in vocab only.
    | cross_check_genes.run(
      fromState: [
          "input": "input",
          "modality": "modality",
          "vocab_file": "model_vocab",
          "var_gene_names": "input_var_gene_names",
          "output": "output",
          "pad_token": "pad_token"
        ],
      toState: ["input": "output"]
    )
    // Bins the data into a fixed number of bins.
    | binning.run(
        fromState: [
            "input": "input",
            "modality": "modality",
            "input_layer": "input_layer",
            "n_input_bins": "n_input_bins",
            "output": "output",
            "seed": "seed"
          ],
        args: [
          "binned_layer": "binned"
        ],
        toState: ["input": "output"]
    )
    // Padding and tokenization of gene count values.
    | pad_tokenize.run(
        fromState: [
          "input": "input",
          "modality": "modality",
          "model_vocab": "model_vocab",
          "input_layer": "binned",
          "var_gene_names": "input_var_gene_names",
          "pad_token": "pad_token",
          "pad_value": "pad_value",
          "max_seq_len": "max_seq_len",
          "output": "output"
        ],
        args: [
          "obsm_gene_tokens": "gene_id_tokens",
          "obsm_tokenized_values": "values_tokenized",
          "obsm_padding_mask": "padding_mask"
        ],
        toState: ["input": "output"]
    )
    // scGPT decoder-based cell type annotation.
    | scgpt_celltype_annotation.run(
        fromState: [
          "model": "model",
          "model_vocab": "model_vocab",
          "model_config": "model_config",
          "label_mapper_key": "label_mapper_key",
          "finetuned_checkpoints_key": "finetuned_checkpoints_key",
          "input": "input",
          "modality": "modality",
          "obs_batch_label": "input_obs_batch_label",
          "pad_token": "pad_token",
          "pad_value": "pad_value",
          "n_input_bins": "n_input_bins",
          "dsbn": "dsbn",
          "batch_size": "batch_size",
          "seed": "seed",
          "output_obs_predictions": "output_obs_predictions",
          "output_obs_probability": "output_obs_probability",
          "output": "workflow_output",
          "output_compression": "output_compression"
        ],
        args: [
          "obsm_gene_tokens": "gene_id_tokens",
          "obsm_tokenized_values": "values_tokenized"
        ],
        toState: {id, output, state -> ["output": output.output]}
    )

  emit:
    output_ch
}
