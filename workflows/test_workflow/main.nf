nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { testProcess } from workflowDir + '/test_workflow/test.nf'

config = readConfig("$workflowDir/test_workflow/config.vsh.yaml")


workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | view { "Input: $it" }
    | map {it -> [it[0], it[1].input]}
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  emit:
  output_ch

  main:
  output_ch = input_ch
    | testProcess

}