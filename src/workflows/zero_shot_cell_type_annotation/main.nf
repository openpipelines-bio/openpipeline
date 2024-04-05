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
            fromState: {id, state ->
            // Annotates the mudata object with highly variable genes.
                [
                "input": state.input,
                // "layer": state.layer,
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
            fromState: {id, state -> 
            [
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
            fromState: {id, state -> 
            [
                "input": state.input,
                "modality": state.modality,
                "model_vocab": state.vocab_file,
                "input_layer": state.binned_layer,
                "gene_name_layer": state.gene_name_layer,
                "pad_token": state.pad_token,
                "pad_value": state.pad_value,
                "max_seq_length": state.max_seq_length,
                "output": state.workflow_output
            ]
            },
            toState: 
            [
                "input": "output",
                "input_obsm_gene_tokens": "output_obsm_gene_tokens",
                "input_obsm_tokenized_values": "output_obsm_padding_mask"
            ]
        )

        | cell_type_inference.run(
            fromState: {id, state -> 
            [
                "input": state.input,
                "modality": state.modality,
                "input_gene_ids": state.input_obsm_gene_tokens,
                "input_values": state.input_obsm_tokenized_values,
                "model": state.model,
                "model_config": state.model_config,
                "model_vocab": state.vocab_file,
                "output": state.workflow_output
            ]
            },
        ) 
        | setState(["output"])

    emit:
        output_ch
}