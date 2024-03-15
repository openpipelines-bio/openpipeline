nextflow.enable.dsl=2

include { split_modalities } from params.rootDir + "/target/nextflow/dataflow/split_modalities/main.nf"
include { split_modalities_test } from params.rootDir + "/target/nextflow/test_workflows/multiomics/split_modalities_test/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
    [
      id: "mouse",
      input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
      publish_dir: "foo/",
      output: "modalities",
      output_types: "types.csv"
    ]
  ])
  | map{ state -> [state.id, state] }
  | split_modalities
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    "Output: $output"
  }
  | map { id, output -> [id, ["input": output.output_types, "mod_dir": output.output]]}
  | split_modalities_test
}