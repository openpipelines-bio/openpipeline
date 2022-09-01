nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { run_singlesample } from "$workflowDir/scrnaseq/singlesample/main.nf"
include { run_multisample } from "$workflowDir/scrnaseq/multisample/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

workflow {
  config = readConfig("$workflowDir/scrnaseq/config.vsh.yaml")
  helpMessage(config)
  viashChannel(params, config)
    | run_singlesample
    | run_multisample
}

workflow singlesample {
  config = readConfig("$workflowDir/scrnaseq/singlesample/config.vsh.yaml")
  helpMessage(config)
  viashChannel(params, config)
    | run_singlesample
}

workflow multisample {
  config = readConfig("$workflowDir/scrnaseq/multisample/config.vsh.yaml")
  helpMessage(config)
  viashChannel(params, config)
    | run_multisample
}
