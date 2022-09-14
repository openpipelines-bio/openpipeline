nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/ingestion/cellranger_mapping/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  main:

  output_ch = input_ch

    // store output value in 3rd slot for later use
    // and transform for concat component
    | map { id, data ->
      new_data = data.clone()
      new_data.remove("output_h5mu")
      new_data.remove("output_raw")
      new_data = new_data + [ output: data.output_raw ]
      
      [id, new_data, data]
    }

    // run count
    | cellranger_count.run(auto: [ publish: true ])

    // split output dir into map
    | cellranger_count_split

    // convert to h5mu
    | map { id, cellranger_outs, orig_data -> 
      [ id, [ input: cellranger_outs.filtered_h5, output: orig_data.output_h5mu ], cellranger_outs ]
    }
    | from_10xh5_to_h5mu.run(auto: [ publish: true ])

    // return output map
    | map { id, h5mu, data -> [ id, data + [h5mu: h5mu] ] }

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
    reference: params.resources_test + "/cellranger_tiny_fastq/cellranger_tiny_ref"
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