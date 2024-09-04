nextflow.enable.dsl=2

include { scgpt_annotation } from params.rootDir + "/target/nextflow/workflows/annotation/scgpt_annotation/main.nf"
include { scgpt_annotation_test } from params.rootDir + "/target/nextflow/test_workflows/annotation/scgpt_annotation_test/main.nf"

workflow test_wf {
    resources_test = file("${params.rootDir}/resources_test/scgpt")

    output_ch = Channel.fromList([
        [
            id: "simple_execution_test",
            input: resources_test.resolve("test_resources/Kim2020_Lung_subset_preprocessed.h5mu"),
            model: resources_test.resolve("finetuned_model/best_model.pt"),
            model_config: resources_test.resolve("source/args.json"),
            model_vocab: resources_test.resolve("source/vocab.json"),
            finetuned_checkpoints_key: "model_state_dict",
            label_mapper_key: "id_to_class",
            input_layer: "log_normalized",
            obs_batch_label: "sample",
            // change default to reduce resource requirements
            n_hvg: 400,
            seed: 1
        ]
    ])
    | map{ state -> [state.id, state] }
    | scgpt_annotation
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test")

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

      "Output: $output"
    }
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert output_list.collect{it[0]} == ["simple_execution_test"]
    }
    | scgpt_annotation_test.run(
        fromState: ["input": "output"]
    )
    | toList()
    | view { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
    }
}
