nextflow.enable.dsl=2

include { zero_shot_cell_type_annotation } from params.rootDir + "/target/nextflow/workflows/zero_shot_cell_type_annotation/main.nf"

workflow test_wf {
    resources_test = file("${params.rootDir}/resources_test")

    output_ch = Channel.fromList([
        [
            // TODO: to be changed to full size/4000 obs dataset
            id: "Kim2020_cpu_run_default",
            input: resources_test.resolve("Kim2020_Lung_cpu_run.h5mu"),
            model: resources_test.resolve("best_model.pt"),
            model_config: resources_test.resolve("args.json"),
            model_vocab: resources_test.resolve("vocab.json"),
            output: "annotation_default.h5mu"
        ],
        [
            // TODO: to be changed to full size/4000 obs dataset
            id: "Kim2020_cpu_run_set_params",
            input: resources_test.resolve("Kim2020_Lung_cpu_run.h5mu"),
            model: resources_test.resolve("best_model.pt"),
            model_config: resources_test.resolve("args.json"),
            model_vocab: resources_test.resolve("vocab.json"),
            gene_name_layer: "gene_name",
            predicted_cell_type_id: "predicted_cell_type",
            pad_token: "<pad>",
            dsbn: True,
            pad_value: -2,
            n_cls: 8,
            n_input_bins: 51,
            batch_size: 64,
            output: "annotation_set_params.h5mu"
        ]
    ])
    | map{ state -> [state.id, state] }
    | zero_shot_cell_type_annotation
    | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
    }
    | toSortedList{a, b -> a[0] <=> b[0]}
    | map { output_list ->
        assert output_list.size() == 2 : "output channel should contain two events"
        println "output_list: $output_list"
        assert output_list.collect{it[0]} == ["Kim2020_cpu_run_default", "Kim2020_cpu_run_set_params"] : "Output ID should be same as input ID"
        assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["mitochondrial_test.final.h5mu", "simple_execution_test.final.h5mu"] as Set)

    }
}