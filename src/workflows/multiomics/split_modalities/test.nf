nextflow.enable.dsl=2

include { split_modalities } from params.rootDir + "/target/nextflow/dataflow/split_modalities/main.nf"
include { split_modalities_test } from params.rootDir + "/target/nextflow/test_workflows/multiomics/split_modalities_test/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  input_ch = Channel.fromList([
    [
      id: "mouse",
      input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
      publish_dir: "foo/",
      output: "modalities",
      output_types: "types.csv"
    ]
  ])
  output_ch = input_ch
  | map { state -> [state.id, state]}
  | split_modalities.run(
    toState: { id, output, state -> output + [orig_input: state.input] }
  )
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    "Output: $output"
  }

  | split_modalities_test.run(
    fromState: {id, state ->
      [ input: state.output_types,
        mod_dir: state.output,
        orig_input: state.orig_input]}
  )
}