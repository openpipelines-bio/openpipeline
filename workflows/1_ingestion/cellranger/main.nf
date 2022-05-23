nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_mkfastq } from targetDir + "/demux/cellranger_mkfastq/main.nf"
include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"
include { from_10xh5_to_h5ad } from targetDir + "/convert/from_10xh5_to_h5ad/main.nf"

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
  --reference      Path to a Cell Ranger reference (required).
  --publishDir     Path to an output directory (required).
  
Parameters (Batch mode):
  --csv            A csv file containing columns 'id', 'input', 'sample_sheet', 'reference' (required).
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
    input_ch = Channel.value( params.subMap(["id", "input", "sample_sheet", "reference"]) )
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
      // process input
      if (li.containsKey("reference") && li.reference) {
        reference_path = file(li.reference)
      } else {
        exit 1, paramExists("csv") ? 
          "ERROR: The provided csv file should contain a 'reference' column" : 
          "ERROR: Please specify an '--reference' parameter"
      }

      // process id
      if (li.containsKey("id") && li.id) {
        id_value = li.id
      } else if (!paramExists("csv")) {
        id_value = "run"
      } else {
        exit 1, "ERROR: The provided csv file should contain an 'id' column"
      }
      [ id_value, [ input: input_path, sample_sheet: sample_sheet_path, reference: reference_path ] ]
    }
    | view { "Input: $it" }
    | run_wf
    | publish.run(
      map: { [ it[0], [ input: it[1], output: it[0] ] ] },
      auto: [ publish: true ]
    )
    | view { "Output: ${params.publishDir}/${it[1]}" }
}

/* Cell Ranger - common workflow
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

  auto = [ publish: paramExists("publishDir"), transcript: paramExists("publishDir") ]
  auto_nopub = [ transcript: paramExists("publishDir") ]

  main:
  output_ch = input_ch

    // run mkfastq
    | map { id, data -> [ id, data.subMap(["input", "sample_sheet"]), data ]}
    | cellranger_mkfastq.run(auto: auto)

    // run count
    | map { id, fastq, data -> [ id, data.subMap("reference") + [ input: fastq ], [ fastq: fastq ] }
    | cellranger_count.run(auto: auto)

    // split output dir into map
    | cellranger_count_split.run(auto: auto_nopub)

    // convert to h5ad
    | map { id, cellranger_outs, data -> [ id, cellranger_outs.filtered_h5, data + cellranger_outs ] }
    | from_10xh5_to_h5ad.run(auto: auto)

    // convert to h5mu
    | map { id, h5ad, data -> [ id, data.filtered_h5, data + [h5ad: h5ad] ] }
    | from_10xh5_to_h5mu.run(auto: auto)

    // return output map
    | map { id, h5mu, data -> [ id, data + [h5mu: h5mu] ] }

  emit:
  output_ch
}


/* Cell Ranger - Integration testing
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