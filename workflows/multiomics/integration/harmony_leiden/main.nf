
nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { harmonypy } from targetDir + '/integrate/harmonypy/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'
include { move_obsm_to_obs } from targetDir + '/metadata/move_obsm_to_obs/main.nf'
include { readConfig; helpMessage; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/harmony_leiden/config.vsh.yaml")

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

    // run harmonypy
    | harmonypy.run(
      fromState: { id, state -> 
        [
          "input": state.input,
          "modality": state.modality,
          "obsm_input": state.embedding,
          "obs_covariates": state.obs_covariates,
          "obsm_output": state.obsm_integrated,
          "theta": state.theta
        ]
      },
      toState: { id, output, state ->
        state + ["input": output]
      }
    )
    
    // run knn
    | find_neighbors.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": state.modality,
          "uns_output": state.uns_neighbors,
          "obsp_distances": state.obsp_neighbor_distances,
          "obsp_connectivities": state.obsp_neighbor_connectivities,
          "obsm_input": state.obsm_integrated
        ]
      },
      toState: { id, output, state ->
        state + ["input": output]
      }
    )

    // run leiden clustering
    | leiden.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": state.modality,
          "obsp_connectivities": state.obsp_neighbor_connectivities,
          "obsm_name": state.obs_cluster,
          "resolution": state.leiden_resolution
        ]
      },
      toState: { id, output, state ->
        state + ["input": output]
      }
    )
    
    // run umap
    | umap.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": state.modality,
          "obsm_input": state.obsm_integrated,
          "obsm_output": state.obsm_umap
        ]
      },
      toState: { id, output, state ->
        state + ["input": output]
      }
    )
    
    // move obsm to obs
    | move_obsm_to_obs.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "modality": state.modality,
          "obsm_key": state.obs_cluster,
          "output": state.output,
          "output_compression": "gzip"
        ]
      },
      // simply replace the state with the output
      toState: { id, output, state ->
        output
      },
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
    param_list: [
      [
        id: "foo",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        layer: "log_normalized",
        obs_covariates: "sample_id",
        embedding: "X_pca",
        leiden_resolution: [1, 0.25],
        output: "foo.final.h5mu"
      ]
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
      assert output_list.size() == 1 : "output channel should contain 1 event"
      assert (output_list.collect({it[0]}) as Set).equals(["foo"] as Set): "Output ID should be same as input ID"
      assert (output_list.collect({it[1].getFileName().toString()}) as Set).equals(["foo.final.h5mu"] as Set)
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}