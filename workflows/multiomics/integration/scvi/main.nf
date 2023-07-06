
nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { scvi } from targetDir + '/integrate/scvi/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'

include { readConfig; helpMessage; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/scvi/config.vsh.yaml")

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
    // split params for downstream components
    | setWorkflowArguments(
      scvi: [
        "input": "input",
        "obs_batch": "obs_batch",
        "obsm_output": "obsm_output",
        "early_stopping": "early_stopping",
        "early_stopping_monitor": "early_stopping_monitor",
        "early_stopping_patience": "early_stopping_patience",
        "early_stopping_min_delta": "early_stopping_min_delta",
        "max_epochs": "max_epochs",
        "reduce_lr_on_plateau": "reduce_lr_on_plateau",
        "lr_factor": "lr_factor",
        "lr_patience": "lr_patience",
        "output_model": "output_model",
        "modality": "modality"
      ],
      neighbors: [
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_output", // use output from scvi as input for neighbors,
        "modality": "modality"
      ],
      umap: [ 
        "uns_neighbors": "uns_neighbors",
        "output": "output",
        "obsm_output": "obsm_umap",
        "modality": "modality",
        "output": "output"
      ]
    )
    | getWorkflowArguments(key: "scvi")
    | scvi
    | pmap {id, arguments, other_arguments -> 
      def input = arguments.output
      def new_arguments = arguments.clone()
      new_arguments.removeAll({k, v -> ["output", "model_output"].contains(k)})
      return [id, new_arguments + ["input": input], other_arguments]
    }
    | getWorkflowArguments(key: "neighbors")
    | find_neighbors
    | getWorkflowArguments(key: "umap")
    | umap.run(
        args: [ output_compression: "gzip" ],     
        auto: [ publish: true ]
    )
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
    | map {list -> list + [test_passthrough: "test"]} 
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