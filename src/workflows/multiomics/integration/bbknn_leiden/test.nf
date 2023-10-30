nextflow.enable.dsl=2

include { bbknn_leiden } from params.rootDir + "/target/nextflow/workflows/multiomics/integration/bbknn_leiden/main.nf"

workflow test_wf {
  resources_test = file("${params.rootDir}/resources_test")

  output_ch =
    Channel.fromList([
      [
        id: "foo",
        input: resources_test.resolve("/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        layer: "log_normalized"
      ]
    ])
    | map{ state -> [state.id, state] }
    | bbknn_leiden
    | view { tup ->
      assert tup.size() == 2 : "outputs should contain two elements; [id, output]"

      // check id
      def id = tup[0]
      assert id == "foo" : "ID should be 'foo'. Found: ${id}"

      // check output
      def output = tup[1]
      assert output instanceof Map: "Output should be a map. Found: ${output}"
      assert "output" in output : "Output should contain key 'output'. Found: ${output}"

      // check h5mu
      def output_h5mu = output.output
      assert output_h5mu.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output}"

      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}

workflow test_wf2 {
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
      [
        id: "foo",
        input: resources_test.resolve("/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        layer: "log_normalized",
        leiden_resolution: []
      ]
    ])
    | map{ state -> [state.id, state] }
    | bbknn_leiden
    | view { tup ->
      assert tup.size() == 2 : "outputs should contain two elements; [id, output]"

      // check id
      def id = tup[0]
      assert id == "foo" : "ID should be 'foo'. Found: ${id}"

      // check output
      def output = tup[1]
      assert output instanceof Map: "Output should be a map. Found: ${output}"
      assert "output" in output : "Output should contain key 'output'. Found: ${output}"

      // check h5mu
      def output_h5mu = output.output
      assert output_h5mu.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output}"

      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}