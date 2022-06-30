nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { pca } from targetDir + '/dimred/pca/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { leiden } from targetDir + '/cluster/leiden/main.nf'

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/viash_workflow_helper.nf"

config = readConfig("$workflowDir/integration/multimodal_integration/config.vsh.yaml")

workflow {
  params.testing = false

  helpMessage(params, config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | pca
    | find_neighbors
    | leiden
    | umap.run(
      auto: [ publish: ! params.testing ]
    )

  emit:
  output_ch
}

workflow test_wf {
  params.testing = true
  
  output_ch =
    Channel.value(
      [
        "foo",
        file(params.rootDir + "/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_ums.h5mu"),
        params
      ]
    )
    | view { "Input: [${it[0]}, ${it[1]}, params]" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain three elements; [id, file, params]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[1]}"
      "Output: [${output[0]}, ${output[1]}, params]"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}