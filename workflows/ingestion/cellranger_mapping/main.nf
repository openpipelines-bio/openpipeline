nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"

include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/ingestion/cellranger_mapping/config.vsh.yaml")

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
      cellranger_count: [
        "input": "input",
        "expect_cells": "expect_cells",
        "chemistry": "chemistry",
        "secondary_analysis": "secondary_analysis",
        "generate_bam": "generate_bam",
        "include_introns": "include_introns"
      ],
      from_10xh5_to_h5mu: [ 
        "output": "output_h5mu",
        "obsm_metrics": "obsm_metrics",
        "output_type": "output_type",
      ]
    )

    | getWorkflowArguments(key: "cellranger_count")
    | cellranger_count.run(auto: [ publish: true ])

    // split output dir into map
    | cellranger_count_split
    // convert to h5mu
    | pmap { id, output_data, other_args -> 
      input_data = other_args.from_10xh5_to_h5mu.output_type == "filtered" ? 
        output_data.filtered_h5 : output_data.raw_h5
      // combine new data for from_10xh5_to_h5mu
      new_data =
        [ 
          input: input_data, 
          input_metrics_summary: output_data.metrics_summary
        ] +
        other_args.from_10xh5_to_h5mu

      // store output to fourth field to return as output
      [ id, new_data, other_args, output_data ]
    }
    | from_10xh5_to_h5mu.run(
      auto: [ publish: true ],
      args: [ output_compression: "gzip" ]
    )
    
    // return output map
    | pmap { id, data, other_args, output_data ->
      [ id, output_data + [h5mu: data] ]
    }

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/cellranger_tiny_fastq/cellranger_tiny_fastq",
    reference: params.resources_test + "/cellranger_tiny_fastq/cellranger_tiny_ref",
    output_type: "filtered",
  ]

  output_ch =
    channelFromParams(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}