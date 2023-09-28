nextflow.enable.dsl=2

targetDir = params.rootDir + "/target/nextflow"

include { demux; readConfig; channelFromParams } from targetDir + "/workflows/ingestion/demux/main.nf"

thisConfig = readConfig(moduleDir.resolve("config.vsh.yaml"))

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
    channelFromParams(params, thisConfig)
    | demux
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, state]"

      def id = output[0]
      assert id.contains("_test")

      def state = output[1]
      assert state.containsKey("output_fastq") : "State should contain output_fastq"
      assert state.output_fastq.isDirectory() : "output_fastq should be a directory."
      assert state.containsKey("output_fastqc") : "State should contain output_fastqc"
      assert state.output_fastqc.isDirectory() : "output_fastqc should be a directory."
      assert state.containsKey("output_multiqc") : "State should contain output_multiqc"
      assert state.output_multiqc.isDirectory() : "output_multiqc should be a directory."

      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 3 : "There should be three outputs"
    }
}
