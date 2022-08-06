nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { bcl_convert } from targetDir + "/demux/bcl_convert/main.nf"
include { bcl2fastq } from targetDir + "/demux/bcl2fastq/main.nf"
include { fastqc } from targetDir + "/qc/fastqc/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$projectDir/config.vsh.yaml")

params.demultiplexer = "mkfastq"

workflow {
  helpMessage(params, config)

  demultiplexers = [ "mkfastq", "bclconvert" ]

  if (!demultiplexers.contains(params.demultiplexer)) {
    exit 1, "ERROR: Please provide a demultiplexer that is one of ${demultiplexers}"
  }

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  mkfastq_ch = input_ch
    | filter{ (it[1].demultiplexer ?: params.demultiplexer) == "mkfastq" }
    | cellranger_mkfastq.run(
        [
          directives: [ publishDir: "${params.publishDir}/fastq" ]
        ]
      )

  bcl_convert_ch = input_ch
    | filter{ (it[1].demultiplexer ?: params.demultiplexer)  == "bclconvert" }
    | bcl_convert.run(
        [
          directives: [ publishDir: "${params.publishDir}/fastq" ]
        ]
      )

  bcl2fastq_ch = input_ch
    | filter{ (it[1].demultiplexer ?: params.demultiplexer)  == "bcl2fastq" }
    | bcl2fastq.run(
        [
          directives: [ publishDir: "${params.publishDir}/fastq" ]
        ]
      )

  all_ch = mkfastq_ch | mix( bcl_convert_ch ) | mix ( bcl2fastq_ch )

  all_ch
      /* | map{ [ it[0], [ mode: "dir", input: it[1] ] ] } */
      | fastqc.run(
          [
            mapData: { it -> [ mode: "dir", input: it ] },
            directives: [ publishDir: "${params.publishDir}/qc" ]
          ]
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
      cores: 2,
      memory: 5,
      demultiplexer: "mkfastq"
    ],
    [
      id: "bcl-convert_test",
      input: params.resources_test + "/demultiplication/bcl-updt",
      sample_sheet: params.resources_test + "/demultiplication/bcl-updt/SampleSheet.csv",
      cores: 2,
      memory: 5,
      demultiplexer: "bclconvert"
    ],
    [
      id: "bcl2fastq_test",
      input: params.resources_test + "/cellranger_tiny_bcl/bcl",
      sample_sheet: params.resources_test + "/cellranger_tiny_bcl/bcl/sample_sheet.csv",
      cores: 2,
      memory: 5,
      demultiplexer: "bcl2fastq",
      ignore_missing: true
    ]
  ]

  output_ch =
    viashChannel(params, config)
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
