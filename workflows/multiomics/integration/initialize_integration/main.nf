nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { pca } from targetDir + '/dimred/pca/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'

include { readConfig; helpMessage; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/initialize_integration/config.vsh.yaml")

workflow initialize_integration {
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
  output_ch = input_ch
    | preprocessInputs("config": config)
    // split params for downstream components
    | setWorkflowArguments(
      pca: [
        "input": "input", 
        "obsm_output": "obsm_pca",
        "var_input": "var_pca_feature_selection",
        "modality": "modality",
        "overwrite": "pca_overwrite"
      ],
      neighbors: [
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "modality": "modality"
      ],
      umap: [ 
        "uns_neighbors": "uns_neighbors",
        "output": "output",
        "obsm_output": "obsm_umap",
        "modality": "modality"
      ]
    )
    | getWorkflowArguments(key: "pca")
    | pca
    | pmap {id, data ->
      def new_data = ["input": data.output]
      [id, new_data]
    }
    | getWorkflowArguments(key: "neighbors")
    | find_neighbors
    | pmap {id, data ->
      def new_data = ["input": data.output]
      [id, new_data]
    }
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
        output: "foo.final.h5mu"
      ],
      [
        id: "foo2",
        input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
        layer: "",
        output: "foo2.final.h5mu"
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain two events"
      assert (output_list.collect({it[0]}) as Set).equals(["foo", "foo2"] as Set): "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["foo.final.h5mu", "foo2.final.h5mu" ] as Set)
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}