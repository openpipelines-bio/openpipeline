nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { clr } from targetDir + '/transform/clr/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'

include { readConfig; helpMessage; channelFromParams; preprocessInputs } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/prot_multisample/config.vsh.yaml")

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
    // add the id to the arguments
    | pmap { id, data ->
      def new_data = data + [ input_id: data.sample_id ]
      [id, new_data]
    }
    | setWorkflowArguments (
      "clr": [:],
      "concat": [:]
    )
    | getWorkflowArguments(key: "concat")
    | concat.run([ modality: "prot" ])
    | getWorkflowArguments(key: "clr")
    | clr
    | map {list -> [list[0], list[1]] + list.drop(3)}
  emit:
  output_ch
}


workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "adt_samples",
    sample_id: "pbmc",
    input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
  ]

  output_ch =
    channelFromParams(testParams, config)
    | map {list -> list + [test_passthrough: "test"]}
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain three elements; [id, file, passthrough]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "adt_samples" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}