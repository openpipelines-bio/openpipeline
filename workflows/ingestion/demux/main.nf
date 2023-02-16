nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { bcl_convert } from targetDir + "/demux/bcl_convert/main.nf"
include { bcl2fastq } from targetDir + "/demux/bcl2fastq/main.nf"
include { fastqc } from targetDir + "/qc/fastqc/main.nf"
include { multiqc } from targetDir + "/qc/multiqc/main.nf"

include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/ingestion/demux/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  commonOptions = [
    args: [ output: "fastq/\$id" ],
    auto: [ publish: true ]
  ]
  preprocessed_ch = input_ch
    | preprocessInputs("config": config)

  mkfastq_ch = preprocessed_ch
    | filter{ it[1].demultiplexer == "mkfastq" }
    | cellranger_mkfastq.run(commonOptions)
    
  bcl_convert_ch = preprocessed_ch
    | filter{ it[1].demultiplexer  == "bclconvert" }
    | bcl_convert.run(commonOptions)

  bcl2fastq_ch = preprocessed_ch
    | filter{ it[1].demultiplexer  == "bcl2fastq" }
    | bcl2fastq.run(commonOptions)

  /* Combine the different demultiplexer channels */
  all_ch =
    mkfastq_ch
      | mix( bcl_convert_ch, bcl2fastq_ch )
      | map {  tup ->
        [tup[0], tup[1].output]
      }

  /* Generate fastqc reports for every sample */
  all_ch
    | fastqc.run(
        [
          args: [ mode: "dir", output: "fastqc/\$id" ],
          auto: [ publish: true ]
        ]
      )

  /* Generate multiqc report */
  all_ch
    | map{ it[1] }
    | toSortedList
    | map{ [ "multiqc", it ] }
    | multiqc.run(
        args: [ output: "multiqc/report" ],
        auto: [ publish: true ]
      )

  output_ch = all_ch

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  params.param_list = [
    [
      id: "mkfastq_test",
      input: params.resources_test + "/cellranger_tiny_bcl/bcl",
      sample_sheet: params.resources_test + "/cellranger_tiny_bcl/bcl/sample_sheet.csv",
      demultiplexer: "mkfastq"
    ],
    [
      id: "bclconvert_test",
      input: params.resources_test + "/cellranger_tiny_bcl/bcl2/",
      sample_sheet: params.resources_test + "/cellranger_tiny_bcl/bcl2/sample_sheet.csv",
      demultiplexer: "bclconvert"
    ],
    [
      id: "bcl2fastq_test",
      input: params.resources_test + "/cellranger_tiny_bcl/bcl",
      sample_sheet: params.resources_test + "/cellranger_tiny_bcl/bcl/sample_sheet.csv",
      demultiplexer: "bcl2fastq",
      ignore_missing: true
    ]
  ]

  output_ch =
    channelFromParams(params, config)
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
      assert output_list.size() == 3 : "There should be three outputs"
    }
}
