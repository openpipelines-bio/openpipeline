nextflow.enable.dsl=2

include { prot_multisample } from params.rootDir + "/target/nextflow/workflows/prot/prot_multisample/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
      [
        id: "adt_samples_axis_0",
        sample_id: "pbmc",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        clr_axis: 0
      ],
      [
        id: "adt_samples_axis_1",
        sample_id: "pbmc",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        clr_axis: 1
      ]
    ])
    | map{ state -> [state.id, state] }
    | prot_multisample
    | view { output ->
      assert output.size() == 2 : "outputs should contain three elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      print "output_list: $output_list"
      assert output_list.size() == 2 : "output channel should contain two events"
      assert output_list.collect({it[0]}).sort() == ["adt_samples_axis_0", "adt_samples_axis_1"] : "Output IDs should be [adt_samples_axis_0, adt_samples_axis_1]"
    }

}
