nextflow.enable.dsl=2

params.rootDir = params.rootDir ?: projectDir + "/../../../.."

include { leiden } from params.rootDir + "/target/_private/nextflow/wrappers/tools/leiden/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "cpu_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        device_type: "cpu",
        resolution: [1.0],
        output_compression: "gzip"
      ]
    ])
    | map { state -> [state.id, state] }
    | leiden
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"
      def id = output[0]
      assert id == "cpu_execution_test" : "Unexpected id: ${id}"
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
      "Output: $output"
    }
    | toSortedList({ a, b -> a[0] <=> b[0] })
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert output_list.collect{ it[0] } == ["cpu_execution_test"]
    }
}

// Exercises the GPU (rapids-singlecell) implementation selected with
// '--device_type gpu'. Requires a CUDA-capable NVIDIA GPU, so this entrypoint is
// listed under 'gpu_tests' in _viash.yaml and skipped by the GitHub Actions CI.
workflow test_gpu_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "gpu_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        device_type: "gpu",
        resolution: [1.0],
        output_compression: "gzip"
      ]
    ])
    | map { state -> [state.id, state] }
    | leiden
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"
      def id = output[0]
      assert id == "gpu_execution_test" : "Unexpected id: ${id}"
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
      "Output: $output"
    }
    | toSortedList({ a, b -> a[0] <=> b[0] })
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert output_list.collect{ it[0] } == ["gpu_execution_test"]
    }
}
