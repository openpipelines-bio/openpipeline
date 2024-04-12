nextflow.enable.dsl=2

include { scgpt_integration } from params.rootDir + "/target/nextflow/workflows/integration/scgpt_integration/main.nf"

workflow test_wf {
    resources_test = file("${params.rootDir}/resources_test/scgpt")

    output_ch = Channel.fromList([
        [
            // TODO: to be changed to full size/4000 obs dataset
            id: "Kim2020_Lung_default",
            input: resources_test.resolve("test_resources/Kim2020_Lung_subset.h5mu"),
            model: resources_test.resolve("source/best_model.pt"),
            model_config: resources_test.resolve("source/args.json"),
            model_vocab: resources_test.resolve("source/vocab.json"),
            output: "annotation_default.h5mu"
        ],
        [
            // TODO: to be changed to full size/4000 obs dataset
            id: "Kim2020_Lung_set_params",
            input: resources_test.resolve("test_resources/Kim2020_Lung_subset.h5mu"),
            model: resources_test.resolve("source/best_model.pt"),
            model_config: resources_test.resolve("source/args.json"),
            model_vocab: resources_test.resolve("source/vocab.json"),
            gene_name_layer: "gene_name",
            input_obs_batch_label: "sample",
            predicted_cell_type_id: "predicted_cell_type",
            pad_token: "<pad>",
            DSBN: "True",
            pad_value: -2,
            n_cls: 8,
            n_input_bins: 51,
            batch_size: 64,
            output: "annotation_set_params.h5mu"
        ]
    ])
    | map{ state -> [state.id, state] }
    | scgpt_integration
    | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
    }
    | toSortedList{a, b -> a[0] <=> b[0]}
    | map { output_list ->
        assert output_list.size() == 2 : "output channel should contain two events"
        println "output_list: $output_list"
        assert output_list.collect{it[0]} == ["Kim2020_Lung_default", "Kim2020_Lung_set_params"] : "Output ID should be same as input ID"
        assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["annotation_default.h5mu", "annotation_set_params.h5mu"] as Set)

    }
}