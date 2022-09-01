nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"
include { filter_with_scrublet } from targetDir + "/filter/filter_with_scrublet/main.nf"
include { do_filter } from targetDir + "/filter/do_filter/main.nf"

include { normalize_total } from targetDir + '/transform/normalize_total/main.nf'
include { log1p } from targetDir + '/transform/log1p/main.nf'
include { filter_with_hvg } from targetDir + '/filter/filter_with_hvg/main.nf'

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"


workflow singlesample {
  config = readConfig("$workflowDir/scrnaseq/singlesample.vsh.yaml")
  helpMessage(config)
  viashChannel(params, config)
    | singlesample
}

workflow multisample {
  config = readConfig("$workflowDir/scrnaseq/multisample.vsh.yaml")
  helpMessage(config)
  viashChannel(params, config)
    | multisample
}

workflow run_singlesample {
  take:
  input_ch

  main:
  output_ch = input_ch
    // cell filtering
    | filter_with_counts
    | do_filter.run(
      args: [ obs_filter: "filter_with_counts" ]
    )
    // doublet calling
    | filter_with_scrublet.run(
      auto: [ publish: true ]
    )
    // TODO: ambient rna correction

  emit:
  output_ch
}

workflow run_multisample {
  take:
  input_ch

  main:
  output_ch = input_ch
    // normalisation
    | normalize_total
    | log1p
    // feature annotation
    | filter_with_hvg.run(
      auto: [ publish: true ]
    )

  emit:
  output_ch
}