nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { pca } from targetDir + '/dimred/pca/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { harmonypy } from targetDir + '/integrate/harmonypy/main.nf'

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { splitParams; combineParams } from workflowDir + "/utils/DataFlowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/config.vsh.yaml")

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
  
    // split params for downstream components
    | splitParams(
      pca: [
        "input": "input", 
        "obsm_output": "obsm_pca"
      ],
      integration: [
        "obsm_input": "obsm_pca",
        "obs_covariates": "obs_covariates",
        "obsm_output": "obsm_integrated"
      ],
      neighbors: [
        "obsm_input": "obsm_integrated", 
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities"
      ],
      clustering: [
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obs_name": "obs_cluster"
      ],
      umap: [ 
        "init_pos": "obsm_integrated",
        "uns_neighbors": "uns_neighbors",
        "output": "output",
        "obsm_output": "obsm_umap"
      ]
    )

    | combineParams("pca")
    | pca

    | combineParams("integration")
    | harmonypy

    | combineParams("neighbors")
    | find_neighbors

    | combineParams("clustering")
    | leiden

    | combineParams("umap")
    | umap.run(
      auto: [ publish: true ]
    )

    // remove splitArgs
    | map { tup ->
      [ tup[0], tup[1] ] + tup.drop(3)
    }

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
    layer: "",
    obs_covariates: "sample_id"
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