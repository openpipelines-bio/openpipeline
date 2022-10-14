nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

// include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { add_id } from targetDir + "/metadata/add_id/main.nf"
include { join_csv } from targetDir + "/metadata/join_csv/main.nf"
include { join_uns_to_obs } from targetDir + "/metadata/join_uns_to_obs/main.nf"
include { concat } from targetDir + "/dataflow/concat/main.nf"
include { from_h5mu_to_h5ad } from targetDir + "/convert/from_h5mu_to_h5ad/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataFlowHelper.nf"

config = readConfig("$workflowDir/custom/ts_concat/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)

    // rename .obs_names and add .obs["sample_id"]
    | pmap{ id, data ->
      def new_data = [
        input_id: id, 
        input: data.input, 
        obs_output: 'sample_id', 
        make_observation_keys_unique: true
      ]
      [ id, new_data, data ]
    }
    | add_id

    // join metadata csv to .obs
    | pmap{ id, file, orig_data ->
      def new_data = [ 
        input: file, 
        input_csv: orig_data.input_metadata,
        obs_key: 'sample_id'
      ]
      [ id, new_data, orig_data ]
    }
    | join_csv

    // temporary step to fix .obs
    | pmap{ id, file, orig_data -> 
      def new_data = [ 
        input: file,
        uns_key: "metrics_cellranger"
      ]
      [ id, new_data, orig_data ]
    }
    | join_uns_to_obs

    // combine into one channel event
    | toSortedList{ a, b -> b[0] <=> a[0] }
    | map { tups -> 
      def new_data = [ 
        input_id: tups.collect{it[0]}, 
        input: tups.collect{it[1]}
      ]
      [ "combined", new_data, tups[0][2] ]
    }

    // concatenate into one h5mu
    | pmap{ id, data, other ->
      [ id, data + [ output: other.output ], other]
    }
    | concat.run(
      auto: [ publish: true ]
    )

    // convert to h5ad
    | pmap{ id, file, other ->
      def output_h5ad = other.output.replaceAll('.h5mu$', ".h5ad")
      [ id, [ input: file, output: output_h5ad ], other]
    }
    | from_h5mu_to_h5ad.run(
      auto: [ publish: true ]
    )
}
