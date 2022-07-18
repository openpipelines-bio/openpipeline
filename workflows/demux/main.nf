nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

println(params.rootDir)

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { bcl_convert } from targetDir + "/demux/bcl_convert/main.nf"
include { publish } from targetDir + "/transfer/publish/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/viash_workflow_helper.nf"

config = readConfig("$projectDir/workflow.yaml")


workflow {
  params.testing = false

  helpMessage(params, config)

  demultiplexers = [ "mkfastq", "bclconvert" ]

  if (!demultiplexers.contains(params.demultiplexer)) {
    exit 1, "ERROR: Please provide a demultiplexer that is one of ${demultiplexers}"
  }

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:

  mkfastq_ch = input_ch
    | filter{ params.demultiplexer == "mkfastq" }
    | cellranger_mkfastq

  bcl_convert_ch = input_ch
    | filter{ params.demultiplexer == "bclconvert" }
    | bcl_convert

  output_ch = mkfastq_ch | mix( bcl_convert_ch )

  emit:
  output_ch
}



    /* | run_wf */
    /* | publish.run( */
    /*   map: { [ it[0], [ input: it[1], output: it[0] ] ] }, */
    /*   auto: [ publish: true ] */
    /* ) */
    /* | view { "Output: ${params.publishDir}/${it[1]}" } */
/* } */


workflow test_wf {
  Channel.value(
      [
        "foo",
        [
          input: file(params.rootDir + "/resources_test/cellranger_tiny_bcl/bcl"),
          sample_sheet: file(params.rootDir + "/resources_test/cellranger_tiny_bcl/bcl/sample_sheet.csv"),
        ],
        params
      ]
    )
    | view { "Input: [${it[0]}, ${it[1]}, params]" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain three elements; [id, file, params]"
      assert output[1].isDirectory() : "Output path should be a directory."
      // todo: check whether output dir contains fastq files
      "Output: [${output[0]}, ${output[1]}, params]"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
}
