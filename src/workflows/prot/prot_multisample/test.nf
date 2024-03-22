nextflow.enable.dsl=2

include { prot_multisample } from params.rootDir + "/target/nextflow/workflows/prot/prot_multisample/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
      [
        id: "adt_samples",
        sample_id: "pbmc",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        axis: 1
      ]
    ])
    | map{ state -> [state.id, state] }
    | prot_multisample
    | view { output ->
      assert output.size() == 2 : "outputs should contain three elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "adt_samples" : "Output ID should be same as input ID"
    }
}
