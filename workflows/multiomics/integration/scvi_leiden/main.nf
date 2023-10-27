
nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { move_obsm_to_obs } from targetDir + '/metadata/move_obsm_to_obs/main.nf'
include { scvi } from targetDir + '/integrate/scvi/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'

include { readConfig; helpMessage; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/scvi_leiden/config.vsh.yaml")

workflow scvi_leiden {
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
  neighbors_ch = input_ch
    | preprocessInputs("config": config)
    | scvi.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "obs_batch": state.obs_batch,
          "obsm_output": state.obsm_output,
          "var_input": state.var_input,
          "early_stopping": state.early_stopping,
          "early_stopping_monitor": state.early_stopping_monitor,
          "early_stopping_patience": state.early_stopping_patience,
          "early_stopping_min_delta": state.early_stopping_min_delta,
          "max_epochs": state.max_epochs,
          "reduce_lr_on_plateau": state.reduce_lr_on_plateau,
          "lr_factor": state.lr_factor,
          "lr_patience": state.lr_patience,
          "output_model": state.output_model,
          "modality": state.modality
       ]
      },
    // use map when viash 0.7.6 is released
    // related to https://github.com/viash-io/viash/pull/515
    //   fromState: [
    //     "input": "input",
    //     "obs_batch": "obs_batch",
    //     "obsm_output": "obsm_output",
    //     "var_input": "var_input",
    //     "early_stopping": "early_stopping",
    //     "early_stopping_monitor": "early_stopping_monitor",
    //     "early_stopping_patience": "early_stopping_patience",
    //     "early_stopping_min_delta": "early_stopping_min_delta",
    //     "max_epochs": "max_epochs",
    //     "reduce_lr_on_plateau": "reduce_lr_on_plateau",
    //     "lr_factor": "lr_factor",
    //     "lr_patience": "lr_patience",
    //     "output_model": "output_model",
    //     "modality": "modality"
    //   ],
      toState: ["input": "output"]
    )
    | find_neighbors.run(
      fromState: [
        "input": "input",
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_output", // use output from scvi as input for neighbors,
        "modality": "modality"
      ],
      toState: ["input": "output"]
    )

  with_leiden_ch = neighbors_ch
    | filter{id, state -> state.leiden_resolution}
    | leiden.run(
      fromState: [
        "input": "input",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution",
        "modality": "modality",
      ],
      toState: ["input": "output"]
    )
    | move_obsm_to_obs.run(
      fromState: [
          "input": "input",
          "obsm_key": "obs_cluster",
          "modality": "modality",
      ],
      toState: ["input": "output"]
    )


  without_leiden_ch = neighbors_ch
    | filter{id, state -> !state.leiden_resolution}

  output_ch = with_leiden_ch.mix(without_leiden_ch)
    | umap.run(
      fromState: {id, state -> [
        "input": state.input,
        "uns_neighbors": state.uns_neighbors,
        "obsm_output": state.obsm_umap,
        "modality": state.modality,
        "output": state.output,
        "output_compression": "gzip"
        ]
      },
      auto: [ publish: true ],
      toState: { id, output, state ->
        [ output: output.output ]
      }
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
        layer: "log_normalized",
        obs_batch: "sample_id",
        max_epochs: 1,
        output: "foo.final.h5mu",
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    // add a test for passthrough
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert (output_list.collect({it[0]}) as Set).equals(["foo"] as Set): "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["foo.final.h5mu"] as Set)
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}


workflow test_wf2 {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
      [
        id: "foo",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        layer: "log_normalized",
        obs_batch: "sample_id",
        max_epochs: 1,
        output: "foo.final.h5mu",
        leiden_resolution: []
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    // add a test for passthrough
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert (output_list.collect({it[0]}) as Set).equals(["foo"] as Set): "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["foo.final.h5mu"] as Set)
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}