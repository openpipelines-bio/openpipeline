nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { bbknn } from targetDir + '/neighbors/bbknn/main.nf'
include { move_obsm_to_obs } from targetDir + '/metadata/move_obsm_to_obs/main.nf'

include { readConfig; helpMessage; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/bbknn_leiden/config.vsh.yaml")

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
  output_ch = input_ch
    | preprocessInputs("config": config)

    // compute bbknn graph
    | bbknn.run(
      fromState: { id, state ->
        [
          input: state.input,
          modality: state.modality,
          obsm_input: state.obsm_input,
          obs_batch: state.obs_batch,
          uns_output: state.uns_output,
          obsp_distances: state.obsp_distances,
          obsp_connectivities: state.obsp_connectivities,
          n_neighbors_within_batch: state.n_neighbors_within_batch,
          n_pcs: state.n_pcs,
          n_trim: state.n_trim
        ]
      },
      // use map when viash 0.7.6 is released
      // related to https://github.com/viash-io/viash/pull/515
      // fromState: [
      //   "input": "input",
      //   "obsm_input": "obsm_input",
      //   "obs_batch": "obs_batch",
      //   "modality": "modality",
      //   "uns_output": "uns_output",
      //   "obsp_distances": "obsp_distances",
      //   "obsp_connectivities": "obsp_connectivities",
      //   "n_neighbors_within_batch": "n_neighbors_within_batch",
      //   "n_pcs": "n_pcs",
      //   "n_trim": "n_trim"
      // ],
      toState: [
        "input": "output"
      ]
    )

    // run leiden on the bbknn graph
    | leiden.run(
      fromState: [
        "input": "input",
        "obsp_connectivities": "obsp_connectivities",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution",
        "modality": "modality"
      ],
      toState: [
        "input": "output"
      ]
    )

    // run umap on the bbknn graph
    | umap.run(
      fromState: [
        "input": "input",
        "uns_neighbors": "uns_output",
        "obsm_output": "obsm_umap",
        "modality": "modality"
      ],
      toState: [
        "input": "output"
      ]
    )

    // move obsm leiden cluster dataframe to obs
    | move_obsm_to_obs.run(
      fromState: { id, state ->
        [
          input: state.input,
          obsm_key: state.obs_cluster,
          modality: state.modality,
          output: state.output,
          output_compression: "gzip"
        ]
      },
      toState: { id, output, state -> 
        [ output: output.output ]
      },
      auto: [publish: true]
    )

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
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        layer: "log_normalized"
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { tup ->
      assert tup.size() == 2 : "outputs should contain two elements; [id, output]"

      // check id
      def id = tup[0]
      assert id == "foo" : "ID should be 'foo'. Found: ${id}"

      // check output
      def output = tup[1]
      assert output instanceof Map: "Output should be a map. Found: ${output}"
      assert "output" in output : "Output should contain key 'output'. Found: ${output}"

      // check h5mu
      def output_h5mu = output.output
      assert output_h5mu.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output}"

      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}