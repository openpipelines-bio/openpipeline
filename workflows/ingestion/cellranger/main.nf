nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/viash_workflow_helper.nf"

config = readConfig("$projectDir/workflow.yaml")

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
  auto_nopub = [ publish: false, transcript: ! params.testing ]

  output_ch = input_ch

    // run mkfastq
    | map { id, data -> [ id, data.subMap(["input", "sample_sheet"]), data ]}
    | cellranger_mkfastq.run(auto: auto)

    // run count
    | map { id, fastq, data -> [ id, data.subMap("reference") + [ input: fastq ], [ fastq: fastq ] ] }
    | cellranger_count.run(auto: auto)

    // split output dir into map
    | cellranger_count_split.run(auto: auto_nopub)

    // convert to h5mu
    | map { id, cellranger_outs, data -> [ id, cellranger_outs.filtered_h5, data + cellranger_outs ] }
    | from_10xh5_to_h5mu.run(auto: auto)

    // return output map
    | map { id, h5mu, data -> [ id, data + [h5mu: h5mu] ] }

  emit:
  output_ch
}


/* Cell Ranger - Integration testing
 */
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
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}