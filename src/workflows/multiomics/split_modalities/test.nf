nextflow.enable.dsl=2

include { split_modalities } from params.rootDir + "/target/nextflow/workflows/multiomics/split_modalities/main.nf"
include { split_modalities_test } from params.rootDir + "/target/_test/nextflow/test_workflows/multiomics/split_modalities_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
      id: "mouse",
      input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
      publish_dir: "foo/",
      output: "modalities",
      output_types: "types.csv"
    ]
  ])
  | map { state -> [state.id, state]}
  | split_modalities.run(
    toState: { id, output, state -> output + [orig_input: state.input] }
  )
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    "Output: $output"
  }

  | split_modalities_test.run(
    fromState: [ 
      "input": "output_types",
      "mod_dir": "output",
      "orig_input": "orig_input",
    ]
  )
}
