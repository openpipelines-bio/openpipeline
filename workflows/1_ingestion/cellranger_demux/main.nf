nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"

include { publish } from targetDir + "/transfer/publish/main.nf" params(params)
include { getChild; paramExists; assertParamExists } from workflowDir + "/utils/utils.nf" params(params)

workflow {
  if (paramExists("help")) {
    log.info """Cell Ranger Demux - CLI workflow

Use cellranger demux to demultiplex sequencing BCL output to FASTQ.
This workflow can be run on a single input or in batch, see below.

Parameters (Single input mode):
  --id             ID of the sample (optional).
  --input          A BCL directory (required).
  --sample_sheet   Sample sheet (required).
  --publishDir     Path to an output directory (required).
  
Parameters (Batch mode):
  --csv            A csv file containing columns 'id', 'input', 'sample_sheet' (required).
  --publishDir     Path to an output directory (required).
"""
    exit 0
  }

  if (paramExists("input") == paramExists("csv")) {
    exit 1, "ERROR: Please provide either an --input parameter or a --csv parameter"
  }
  
  assertParamExists("publishDir", "where output files will be published")

  if (paramExists("csv")) {
    input_ch = Channel.fromPath(params.csv)
      | splitCsv(header: true, sep: ",")
  } else {
    input_ch = Channel.value( params.subMap(["id", "input", "sample_sheet"]) )
  }

  input_ch
    | map { li ->
      // process input
      if (li.containsKey("input") && li.input) {
        input_path = li.input.split(";").collect { path -> 
          file(paramExists("csv") ? getChild(params.csv, path) : path)
        }.flatten()
      } else {
        exit 1, paramExists("csv") ? 
          "ERROR: The provided csv file should contain an 'input' column" : 
          "ERROR: Please specify an '--input' parameter"
      }
      // process input
      if (li.containsKey("sample_sheet") && li.sample_sheet) {
        sample_sheet_path = file(li.sample_sheet)
      } else {
        exit 1, paramExists("csv") ? 
          "ERROR: The provided csv file should contain a 'sample_sheet' column" : 
          "ERROR: Please specify an '--sample_sheet' parameter"
      }

      // process id
      if (li.containsKey("id") && li.id) {
        id_value = li.id
      } else if (!paramExists("csv")) {
        id_value = "run"
      } else {
        exit 1, "ERROR: The provided csv file should contain an 'id' column"
      }
      [ id_value, [ input: input_path, sample_sheet: sample_sheet_path ] ]
    }
    | view { "Input: $it" }
    | run_wf
    | publish.run(
      map: { [ it[0], [ input: it[1], output: it[0] ] ] },
      auto: [ publish: true ]
    )
    | view { "Output: ${params.publishDir}/${it[1]}" }
}

/* Cell Ranger Demux - common workflow
 * 
 * consumed params:
 *   id:                            a sample id for one or more fastq files
 *   data:
 *     input:                       a BCL directory (required).
 *   output                         a publish dir for the output h5ad files
 * output format:               [ id, fastq, params ]
 *   value id:                      a sample id for one or more fastq files
 *   value fastq_dir:               a directory of fastq files
 *   value params:                  the params object, which may already have sample specific overrides
 */
workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | cellranger_mkfastq

  emit:
  output_ch
}


/* Cell Ranger Demux - Integration testing
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