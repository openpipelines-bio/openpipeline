nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/viash_workflow_helper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

workflow {
  params.testing = false

  helpMessage(params, config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  auto = [ publish: ! params.testing, transcript: ! params.testing ]

  output_ch = input_ch
    | cellranger_mkfastq.run(auto: auto)

  emit:
  output_ch
}

workflow test_wf {
  params.testing = true
  
  Channel.value(
      [
        "foo",
        [
          input: file(params.rootDir + "/resources_test/cellranger_tiny_bcl/bcl"),
          sample_sheet: file(params.rootDir + "/resources_test/cellranger_tiny_bcl/bcl/sample_sheet.csv"),
          cores: 2,
          memory: 5
        ]
      ]
    )
    | view { "Input: $it" }
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
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}