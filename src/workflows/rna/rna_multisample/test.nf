nextflow.enable.dsl=2

include { rna_multisample } from params.rootDir + "/target/nextflow/workflows/rna/rna_multisample/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "simple_execution_test",
        input: resources_test.resolve("concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu"),
        output: "concatenated_file.final.h5mu"
      ]
    ])
    | map{ state -> [state.id, state] }
    | rna_multisample
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
      assert output_list.size() == 1 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["simple_execution_test"]
    }
}

// Exercises the GPU (rapids-singlecell) implementations selected with
// '--device_type gpu'. Requires a CUDA-capable NVIDIA GPU, so this entrypoint is
// listed under 'gpu_tests' in _viash.yaml and skipped by the GitHub Actions CI.
// Scaling is enabled explicitly so the GPU 'scale' variant is covered as well.
workflow test_gpu_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "gpu_execution_test",
        input: resources_test.resolve("concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu"),
        output: "concatenated_file.final.h5mu",
        enable_scaling: true,
        device_type: "gpu"
      ]
    ])
    | map{ state -> [state.id, state] }
    | rna_multisample
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      def id = output[0]
      assert id == "gpu_execution_test"

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
      assert output_list.collect{it[0]} == ["gpu_execution_test"]
    }
}
