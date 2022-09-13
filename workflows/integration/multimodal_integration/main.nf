nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { pca } from targetDir + '/dimred/pca/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { harmonypy } from targetDir + '/integrate/harmonypy/main.nf'

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/integration/multimodal_integration/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // store output value in 3rd slot for later use
    // derive key for pca obsm output
    | map { id, data -> 
      new_data = data + [ obsm_output: data.layer + "_pca" ]
      [ id, new_data, data ]
    }
    | pca

    // set obsm key because pca might not be stored in X_pca
    | map { id, file, orig_data -> 
      new_data = [
        input: file, 
        obsm_input: orig_data.layer + "_pca",
        obsm_output: orig_data.layer + "_pca_integrated",
        obs_covariates: [ "sample_id" ] // temporary override
      ]
      [ id, new_data, orig_data ]
    }
    | harmonypy

    | map { id, file, orig_data -> 
      new_data = [
        input: file, 
        obsm_input: orig_data.layer + "_pca"
      ]
      [ id, new_data, orig_data ]
    }
    | find_neighbors
    | leiden

    // retrieve output value
    | map { id, file, orig_data -> 
      [ id, [ input: file ] + orig_data.subMap("output") ]
    }
    | umap.run(
      auto: [ publish: true ]
    )

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/concat/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
    layer: ""
  ]

  output_ch =
    viashChannel(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}