workflow run_wf {

    take:
        input_ch
    
    main:
        output_ch = input_ch
        | map {id, state -> 
            def new_state = state + ["workflow_output": state.output]
            [id, new_state]
        }

        | highly_variable_features_scanpy.run(
        // Annotates the mudata object with highly variable genes.
            fromState: {id, state ->
                [
                    "input": state.input,
                    "modality": state.modality,
                    "n_top_features": state.n_hvg,
                    "layer": state.input_layer,
                    "var_name_filter": "filter_with_hvg",
                    "flavor": "seurat_v3"
                ]
            },
            toState: ["input": "output"]
        )
        | do_filter.run(
        // do_filter does not need a layer argument because it filters all layers
        // from a modality.
        // filters the mudata object based on the HVG
            fromState: {id, state ->
                [
                    "input": state.input,
                    "modality": state.modality
                    "var_filter": "filter_with_hvg"
                ]
            },
            toState: ["input": "output"]
        )

        | cross_check_genes.run(
        // Check whether the genes are part of the provided vocabulary. Subsets for genes present in vocab only.
            fromState: { id, state ->
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
            fromState: {id, state -> 
            [
                "input": state.input,
                "modality": state.modality,
                "input_layer": state.input_layer,
                "n_input_bins": state.n_input_bins,
                "binned_layer": "binned",
                "seed", state.seed,
                "output": state.output
            ]
            },
            toState: ["input": "output"],
        )

        | pad_tokenize.run(
        // Tokenizes the gene names and pads the tokens.
            fromState: {id, state -> 
            [
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

        | cell_type_inference.run(
        // Infers cell types using the model.
            fromState: {id, state -> 
            [
                "input": state.input,
                "modality": state.modality,
                "obsm_gene_tokens": "gene_id_tokens",
                "obsm_tokenized_values": "values_tokenized",
                "model": state.model,
                "model_config": state.model_config,
                "model_vocab": state.model_vocab,
                "obs_batch_label": state.input_obs_batch_label,
                "obs_predicted_cell_type": state.obs_predicted_cell_type,
                "pad_token": state.pad_token,
                "dsbn": state.dsbn,
                "pad_value": state.pad_value,
                "seed", state.seed,
                "n_cls": state.n_cls,
                "n_input_bins": state.n_input_bins,
                "batch_size": state.batch_size,
                "output": state.workflow_output,
                "output_compression": state.output_compression
            ]
            },
        ) 
        | setState(["output"])

    emit:
        output_ch
}