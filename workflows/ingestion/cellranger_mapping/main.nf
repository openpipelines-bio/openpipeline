nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/ingestion/cellranger_mapping/config.vsh.yaml")

// keep track of whether this is an integration test or not
global_params = [ do_publish: true ]

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  auto = [ publish: global_params.do_publish, transcript: global_params.do_publish ]
  auto_nopub = [ publish: false, transcript: global_params.do_publish ]

  output_ch = input_ch

    // run count
    | cellranger_count.run(auto: auto)

    // split output dir into map
    | cellranger_count_split.run(auto: auto_nopub)

    // convert to h5mu
    | map { id, cellranger_outs -> [ id, cellranger_outs.filtered_h5, cellranger_outs ] }
    | from_10xh5_to_h5mu.run(auto: auto)

    // return output map
    | map { id, h5mu, data -> [ id, data + [h5mu: h5mu] ] }

  emit:
  output_ch
}

workflow test_wf {
  // don't publish output
  global_params.do_publish = false

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/cellranger_tiny_fastq/cellranger_tiny_fastq",
    reference: params.resources_test + "/cellranger_tiny_fastq/cellranger_tiny_ref",
    cores: 2,
    memory: 5
  ]

  output_ch =
    viashChannel(testParams, config)
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