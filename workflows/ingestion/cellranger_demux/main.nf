nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { bcl_convert } from targetDir + "/demux/bcl_convert/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/ingestion/cellranger_demux/config.vsh.yaml")

params.demultiplexer = "mkfastq" // default

workflow {
  helpMessage(config)

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

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/cellranger_tiny_bcl/bcl",
    sample_sheet: params.resources_test + "/cellranger_tiny_bcl/bcl/sample_sheet.csv",
    cores: 2,
    memory: 5
  ]

  output_ch =
    viashChannel(testParams, config)
    | view{ "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].isDirectory() : "Output path should be a directory."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
}
