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

config = readConfig("$workflowDir/custom/ts_process_from10xraw/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)

    // map data to reference
    | pmap{ id, data ->
      [ id, data, data]
    }

    // split output dir into map
    | cellranger_count_split.run(
      args: [ bam: null, bai: null ] // don't copy bam
    )

    // convert to h5mu
    | pmap{ id, data, orig_data -> 
      new_data = [ 
        input: data.raw_h5,
        input_metrics_summary: data.metrics_summary
      ]
      [ id, new_data, orig_data ]
    }
    | from_10xh5_to_h5mu

    // run cellbender
    | cellbender_remove_background.run(
      args: [
        min_counts: 1000,
        layer_output: "cellbender"
      ]
    )

    // filter counts
    | filter_with_counts.run(
      args: [
        layer: "cellbender",
        min_genes: 100, 
        min_counts: 1000, 
        do_subset: true
      ]
    )

    | pmap{ id, file, orig_data -> 
      new_data = [ 
        input: file,
        uns_key: "metrics_cellranger",
        output: orig_data.output_h5mu
      ]
      [ id, new_data ]
    }
    | join_uns_to_obs.run(auto: [ publish: true ])
}
