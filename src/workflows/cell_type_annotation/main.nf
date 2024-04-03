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
                ]
            },
            args: [
                output_gene_ids: "output_gene_ids.pt",
                output_values: "output_values.pt",
                output_padding_mask: "output_padding_mask.pt"
            ],
            toState: [
                "input_gene_ids": "output_gene_ids",
                "input_values": "output_values",
                "input_padding_mask": "output_padding_mask"
            ]
        )

        | annotation.run(
            fromState: {id, state -> [
                "input": state.input,
                "modality": state.modality,
                "input_gene_ids": state.input_gene_ids,
                "input_values": state.input_values,
                "input_padding_mask": state.input_padding_mask,
                "model": state.model,
                "model_config": state.model_config,
                "model_vocab": state.vocab_file,
                "output": state.workflow_output
                ]
            },
            // auto: [ publish: true ]
        ) 
        | setState(["output"])

    emit:
        output_ch
}