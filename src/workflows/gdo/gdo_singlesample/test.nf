nextflow.enable.dsl=2

include { gdo_singlesample } from params.rootDir + "/target/nextflow/workflows/gdo/gdo_singlesample/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
      [
        id: "simple_execution_test",
        input: resources_test.resolve("10x_5k_lung_crispr/SC3_v3_NextGem_DI_CRISPR_A549_5K.h5mu"),
        min_counts: 3,
        max_counts: 10000000,
        min_guides_per_cell: 2,
        max_guides_per_cell: 10000000,
        min_cells_per_guide: 2,
        output: "simple_execution_test.final.h5mu"
      ]
    ])
    | map{ state -> [state.id, state] }
    | gdo_singlesample
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toSortedList{a, b -> a[0] <=> b[0]}
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain two events"
      println "output_list: $output_list"
      assert output_list.collect{it[0]} == ["simple_execution_test"] : "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["simple_execution_test.final.h5mu"] as Set)

    }
}