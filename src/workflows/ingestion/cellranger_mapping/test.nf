nextflow.enable.dsl=2

include { cellranger_mapping } from params.rootDir + "/target/nextflow/workflows/ingestion/cellranger_mapping/main.nf"
include { cellranger_mapping_test } from params.rootDir + "/target/_test/nextflow/test_workflows/ingestion/cellranger_mapping_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [  
        id: "foo",
        input: resources_test.resolve("cellranger_tiny_fastq/cellranger_tiny_fastq"),
        reference: resources_test.resolve("cellranger_tiny_fastq/cellranger_tiny_ref"),
        output_type: "filtered",
      ]
    ])
    | map{ state -> [state.id, state] }
    | cellranger_mapping
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }

    | cellranger_mapping_test.run(
      fromState: ["input": "output_h5mu"]
    )

    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
}
