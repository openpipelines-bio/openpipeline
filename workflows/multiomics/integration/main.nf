nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { pca } from targetDir + '/dimred/pca/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { harmonypy } from targetDir + '/integrate/harmonypy/main.nf'

include { readConfig; helpMessage; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  main:
  pca_ch = input_ch
    | preprocessInputs("config": config)
    // split params for downstream components
    | setWorkflowArguments(
      pca: [
        "input": "input", 
        "obsm_output": "obsm_pca",
        "var_input": "var_pca_feature_selection"
      ],
      integration: [
        "obsm_input": "obsm_pca",
        "obs_covariates": "obs_covariates",
        "obsm_output": "obsm_integrated",
        "theta": "rna_theta"
      ],
      neighbors: [
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities"
      ],
      clustering: [
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obs_name": "obs_cluster",
        "resolution": "leiden_resolution",
      ],
      umap: [ 
        "uns_neighbors": "uns_neighbors",
        "output": "output",
        "obsm_output": "obsm_umap"
      ]
    )
    | pmap { id, args, other_args ->
      // If obs_covariates is not set or empty, harmony will not be run
      // In this case, the layer that find_neighbour uses should not originate from harmonypy but from pca
      def obs_cov = other_args.integration.obs_covariates ?: []
      if (obs_cov.empty) {
        new_neighbors = other_args.neighbors + [obsm_input: other_args.pca.obsm_output]
      } else {
        new_neighbors = other_args.neighbors + [obsm_input: other_args.integration.obsm_output]
      }
      new_other = other_args + [neighbors: new_neighbors]
      [id, args, new_other]
    }
    | getWorkflowArguments(key: "pca")
    | pca
    | getWorkflowArguments(key: "integration")

  without_harmony_ch = pca_ch
    | filter{!it[1].obs_covariates || (it[1].obs_covariates && it[1].obs_covariates.empty)}

  with_harmony_ch = pca_ch
    | filter{it[1].obs_covariates && !it[1].obs_covariates.empty}
    | harmonypy

  output_ch = without_harmony_ch.concat(with_harmony_ch)
    | getWorkflowArguments(key: "neighbors")
    | find_neighbors

    | getWorkflowArguments(key: "clustering")
    | leiden

    | getWorkflowArguments(key: "umap")
    | umap.run(
      auto: [ publish: true ],
      args: [ output_compression: "gzip" ]
    )

    // remove splitArgs
    | map { tup ->
      tup.take(2) + tup.drop(3)
    }

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
      [
        id: "foo",
        input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
        layer: "",
        obs_covariates: "sample_id"
      ],
      [
        id: "foo_without_harmony",
        input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
        layer: "",
        obs_covariates: [],
        leiden_resolution: 2,
      ],
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      println "output_list: $output_list"
      assert output_list.size() == 2 : "output channel should contain two events"
      assert (output_list.collect({it[0]}) as Set).equals(["foo_without_harmony", "foo"] as Set): "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}