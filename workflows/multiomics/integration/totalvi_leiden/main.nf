
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

workflow {
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
  output_ch = integrated_ch 
    | setWorkflowArguments(
      neighbors: [
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_output", // use output from scvi as input for neighbors,
        "query_modality": "modality"
      ],
      clustering: [
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution",
        "query_modality": "modality",
      ],
      umap: [ 
        "uns_neighbors": "uns_neighbors",
        "obsm_output": "obsm_umap",
        "query_modality": "modality",
      ],
      move_obsm_to_obs_leiden: [
        "obsm_key": "obs_cluster",
        "query_modality": "modality"
      ],
      "publish": ["output": "output"]
    )
    | getWorkflowArguments(key: "neighbors")
    | find_neighbors
    | getWorkflowArguments(key: "clustering")
    | leiden
    | getWorkflowArguments(key: "umap")
    | umap
    | getWorkflowArguments(key: "move_obsm_to_obs_leiden")
    | move_obsm_to_obs.run(
      args: [ output_compression: "gzip" ],
      auto: [ publish: true ],
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
    // split params for downstream components
    | setWorkflowArguments(
      totalvi: [
        "input": "input",
        "layer": "input_layer",
        "obs_batch": "obs_batch",
        "query_modality": "modality",
        "query_proteins_modality": "modality_proteins",
        "query_model_path": "query_model_path",
        "obsm_normalized_rna_output": "obsm_rna_output",
        "obsm_normalized_protein_output": "obsm_prot_output",
        "reference_model_path": "reference_model_path",
        "reference_modality": "reference_modality",
        "reference_proteins_modality": "reference_proteins_modality",
        "var_input": "var_input",
        "force_retrain": "force_retrain",
        "weight_decay": "weight_decay",
        "max_epochs": "max_epochs",
        "max_query_epochs": "max_query_epochs"
      ],
      neighbors_leiden_umap_rna: [
        "uns_neighbors": "rna_uns_neighbors",
        "obsp_neighbor_distances": "rna_obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "rna_obsp_neighbor_connectivities",
        "obsm_output": "rna_obsm_output",
        "obs_cluster": "rna_obs_cluster",
        "leiden_resolution": "rna_leiden_resolution",
        "modality": "modality",
        "uns_neighbors": "rna_uns_neighbors",
        "obsm_umap": "obsm_umap",
        "modality": "modality"
      ],
      neighbors_leiden_umap_prot: [
        "uns_neighbors": "prot_uns_neighbors",
        "obsp_neighbor_distances": "prot_obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "prot_obsp_neighbor_connectivities",
        "obsm_output": "prot_obsm_output",
        "obs_cluster": "prot_obs_cluster",
        "leiden_resolution": "prot_leiden_resolution",
        "modality": "modality",
        "uns_neighbors": "prot_uns_neighbors",
        "obsm_umap": "obsm_umap",
        "modality": "modality_proteins"

      ],
    )
    | getWorkflowArguments(key: "totalvi")
    | view {"Before totalVI: $it"}
    | totalvi
    | pmap {id, arguments, other_arguments -> 
      def input = arguments.output
      def new_arguments = arguments.clone()
      new_arguments.removeAll({k, v -> ["output", "model_output"].contains(k)})
      return [id, new_arguments + ["input": input], other_arguments]
    }
    | getWorkflowArguments(key: "neighbors_leiden_umap_rna")
    | neighbors_leiden_umap
    | getWorkflowArguments(key: "neighbors_leiden_umap_prot")
    | neighbors_leiden_umap
    | getWorkflowArguments(key: "move_obsm_to_obs_leiden")

    | pmap {id, arguments, other_arguments ->
      return [id, arguments]
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
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        reference: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
        modality_proteins: "prot",
        reference_proteins_modality: "prot",
        layer: "log_normalized",
        var_input: "filter_with_hvg",
        reference_model_path: "totalvi_reference_model",
        query_model_path: "totalvi_query_model",
        max_epochs: 1,
        max_query_epochs: 1,
        obs_batch: "sample_id",
        output: "foo.final.h5mu"
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    // add a test for passthrough
    | map {list -> list + [test_passthrough: "test"]}
    // Avoid duplicate file names
    | pmap {id, arguments -> 
      copy = arguments.reference.copyTo("reference.h5mu")
      arguments.reference = copy
      [id, arguments]
    }
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain two elements; [id, file, passthrough]"
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