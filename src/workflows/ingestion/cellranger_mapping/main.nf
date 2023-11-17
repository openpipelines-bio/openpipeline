nextflow.enable.dsl=2

workflowDir = params.rootDir + "/src/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"

include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/ingestion/cellranger_mapping/config.vsh.yaml")

workflow cellranger_mapping {
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
    | cellranger_count.run(
      fromState: [
        "input": "input",
        "output": "output_raw",
        "expect_cells": "expect_cells",
        "chemistry": "chemistry",
        "secondary_analysis": "secondary_analysis",
        "generate_bam": "generate_bam",
        "include_introns": "include_introns",
        "reference": "reference"
      ],
      toState: [
        "input": "output",
        "output_raw": "output"
      ],
      auto: [ publish: true ]
    )
    // split output dir into map
    | cellranger_count_split.run(
      fromState: {id, state -> 
        def stateMapping = [
          "input": state.input,
        ]
        outputType = state.output_type == "raw" ? "raw_h5" : "filtered_h5"
        stateMapping += [outputType: "\$id.\$key.${outputType}.h5"]
        stateMapping += ["metrics_summary": "\$id.\$key.metrics_summary.csv"]
        return stateMapping
      },
      toState: {id, output, state -> 
        def outputFile = state.output_type == "raw" ? output.raw_h5 : output.filtered_h5
        def newState = state + [ "input": outputFile ] 
        return newState
      }
    )
    // convert to h5mu
    | from_10xh5_to_h5mu.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "output_compression": "gzip",
          "output": state.output_h5mu,
          "uns_metrics": state.uns_metrics,
          "input_metrics_summary": state.metrics_summary
        ]
      },
      toState: { id, output, state ->
        [
          "output_raw": state.output_raw,
          "output_h5mu": output.output
        ]
      },
      auto: [ publish: true ],
    )

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