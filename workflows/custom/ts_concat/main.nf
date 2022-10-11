nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { cellbender_remove_background } from targetDir + "/correction/cellbender_remove_background/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"
include { add_id } from targetDir + "/metadata/add_id/main.nf"
include { join_csv } from targetDir + "/metadata/join_csv/main.nf"
include { join_uns_to_obs } from targetDir + "/metadata/join_uns_to_obs/main.nf"
include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"
include { concat } from targetDir + "/dataflow/concat/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataFlowHelper.nf"

config = readConfig("$workflowDir/custom/ts_concat/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)

    // rename .obs_names and add .obs["sample_id"]
    | pmap{ id, data ->
      new_data = [
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
      new_data = [ 
        input: file, 
        input_csv: orig_data.input_metadata,
        obs_key: 'sample_id'
      ]
      [ id, new_data ]
    }
    | join_csv

    // combine into one mudata
    | toSortedList{ a, b -> b[0] <=> a[0] }
    | map { tups -> 
      new_data = [ 
        input_id: tups.collect{it[0]}, 
        input: tups.collect{it[1]},
        output: tups[0][2].output
      ]
      [ "combined", new_data ]
    }
    | concat.run(
      auto: [ publish: true ]
    )
}
