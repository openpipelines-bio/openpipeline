workflow run_wf {

    take:
        input_ch
    
    main:
        output_ch = input_ch
        | map {id, state -> 
            def new_state = state + ["workflow_output": state.output]
            [id, new_state]
        }

        | pad_tokenize.run(
            fromState: {id, state -> [
                "input": state.input,
                "modality": state.modality,
                "model_vocab": state.vocab_file,
                "input_layer": state.input_layer,
                "gene_name_layer": state.gene_name_layer,
                "pad_token": state.pad_token,
                "pad_value": state.pad_value,
                "max_seq_length": state.max_seq_length,
                "output": state.workflow_output
                ]
            },
            toState: [
                "input": "output",
                "input_obsm_gene_tokens": "output_obsm_gene_tokens",
                "input_obsm_tokenized_values": "output_obsm_padding_mask",            ]
        )

        | annotation.run(
            fromState: {id, state -> [
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