
nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { move_obsm_to_obs } from targetDir + '/metadata/move_obsm_to_obs/main.nf'
include { totalvi } from targetDir + '/integrate/totalvi/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'
include { publish } from targetDir + '/transfer/publish/main.nf'

include { readConfig; helpMessage; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/totalvi_leiden/config.vsh.yaml")

workflow totalvi_leiden {
  helpMessage(config)

  channelFromParams(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow neighbors_leiden_umap {
  take:
  integrated_ch

  main:
  neighbors_ch = integrated_ch
    | find_neighbors.run(
      fromState: [
        "input": "input",
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_output", // use output from scvi as input for neighbors,
        "query_modality": "modality"
      ],
      toState: ["input": "output"]
    )

  with_leiden_ch = neighbors_ch
    | filter{list -> list[1].leiden_resolution}
    | leiden.run(
      fromState: [
        "input": "input",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution",
        "query_modality": "modality",
      ],
      toState: ["input": "output"]
    )
    | move_obsm_to_obs.run(
      fromState: [
        "input": "input",
        "obsm_key": "obs_cluster",
        "query_modality": "modality",
      ],
      toState: ["input": "output"]
    )

  without_leiden_ch = neighbors_ch
    | filter{list -> !list[1].leiden_resolution}

  output_ch = with_leiden_ch.mix(without_leiden_ch)
    | umap.run(
      fromState: [
          "input": "input",
          "uns_neighbors": "uns_neighbors",
          "obsm_output": "obsm_umap",
          "query_modality": "modality",
        ],
      toState: ["output": "output"]

    )

  emit:
  output_ch
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs("config": config)
    // Avoid conflict with other output arguments
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | totalvi.run(
      fromState: [
        "input": "input",
        "layer": "layer",
        "obs_batch": "obs_batch",
        "query_modality": "modality",
        "query_proteins_modality": "prot_modality",
        "query_model_path": "query_model_path",
        "obsm_normalized_rna_output": "rna_obsm_output",
        "obsm_normalized_protein_output": "prot_obsm_output",
        "reference_model_path": "reference_model_path",
        "reference_modality": "rna_reference_modality",
        "reference_proteins_modality": "prot_reference_modality",
        "var_input": "var_input",
        "force_retrain": "force_retrain",
        "weight_decay": "weight_decay",
        "max_epochs": "max_epochs",
        "max_query_epochs": "max_query_epochs",
        "reference": "reference"
      ],
      toState: ["input": "output"]
    )
    | map { id, state -> // for gene expression
      stateMapping = [
        "input": "input",
        "uns_neighbors": "rna_uns_neighbors",
        "obsp_neighbor_distances": "rna_obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "rna_obsp_neighbor_connectivities",
        "obsm_output": "rna_obsm_output",
        "obs_cluster": "rna_obs_cluster",
        "leiden_resolution": "rna_leiden_resolution",
        "uns_neighbors": "rna_uns_neighbors",
        "obsm_umap": "obsm_umap",
        "modality": "modality"
      ]
      def new_state = stateMapping.collectEntries{newKey, origKey ->
        [newKey, state[origKey]]
      }
      [id, new_state, state]
    }
    | neighbors_leiden_umap
    | map { id, state, orig_state -> // for ADT
      stateMapping = [
        "uns_neighbors": "prot_uns_neighbors",
        "obsp_neighbor_distances": "prot_obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "prot_obsp_neighbor_connectivities",
        "obsm_output": "prot_obsm_output",
        "obs_cluster": "prot_obs_cluster",
        "leiden_resolution": "prot_leiden_resolution",
        "uns_neighbors": "prot_uns_neighbors",
        "obsm_umap": "obsm_umap",
        "modality": "prot_modality",
        "workflow_output": "workflow_output"
      ]
      def new_state = stateMapping.collectEntries{newKey, origKey ->
        [newKey, orig_state[origKey]]
      }
      [id, new_state + ["input": state.output]]
    }
    | neighbors_leiden_umap
    | publish.run(
      fromState: { id, state -> [
          "input": state.output,
          "output": state.workflow_output,
          "compression": "gzip"
        ]
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
        reference: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        prot_modality: "prot",
        prot_reference_modality: "prot",
        var_input: "filter_with_hvg",
        reference_model_path: "totalvi_reference_model",
        query_model_path: "totalvi_query_model",
        max_epochs: 1,
        max_query_epochs: 1,
        output: "foo.final.h5mu"
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    // add a test for passthrough
    // Avoid duplicate file names
    | pmap {id, arguments -> 
      copy = arguments.reference.copyTo("reference.h5mu")
      arguments.reference = copy
      [id, arguments]
    }
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
        reference: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        prot_modality: "prot",
        prot_reference_modality: "prot",
        var_input: "filter_with_hvg",
        reference_model_path: "totalvi_reference_model",
        query_model_path: "totalvi_query_model",
        max_epochs: 1,
        max_query_epochs: 1,
        output: "foo.final.h5mu",
        rna_leiden_resolution: []
      ],
      [
        id: "foo2",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        reference: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        prot_modality: "prot",
        prot_reference_modality: "prot",
        var_input: "filter_with_hvg",
        reference_model_path: "totalvi_reference_model",
        query_model_path: "totalvi_query_model",
        max_epochs: 1,
        max_query_epochs: 1,
        output: "foo2.final.h5mu",
        prot_leiden_resolution: []
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    // add a test for passthrough
    // Avoid duplicate file names
    | pmap {id, arguments -> 
      copy = arguments.reference.copyTo("reference.h5mu")
      arguments.reference = copy
      [id, arguments]
    }
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert (output_list.collect({it[0]}) as Set).equals(["foo", "foo2"] as Set): "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["foo.final.h5mu", "foo2.final.h5mu"] as Set)
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}