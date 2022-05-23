nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellranger_count } from targetDir + "/mapping/cellranger_count/main.nf"
include { cellranger_count_split } from targetDir + "/mapping/cellranger_count_split/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"
include { from_10xh5_to_h5ad } from targetDir + "/convert/from_10xh5_to_h5ad/main.nf"

include { publish } from targetDir + "/transfer/publish/main.nf" params(params)
include { getChild; paramExists; assertParamExists } from workflowDir + "/utils/utils.nf" params(params)

workflow {
  if (paramExists("help")) {
    log.info """Cell Ranger Mapping - CLI workflow

A workflow for running a Cell Ranger Mapping workflow.
This workflow can be run on a single input or in batch, see below.

Parameters (Single input mode):
  --id             ID of the sample (optional).
  --input          One or more fastq paths, separated with semicolons (required).
                   Paths may be globs. Example: path/to/dir/**.fastq
  --reference      Path to a Cell Ranger reference (required).
  --publishDir     Path to an output directory (required).
  
Parameters (Batch mode):
  --csv            A csv file containing columns 'id', 'input' (required).
  --reference      Path to a Cell Ranger reference (required).
  --publishDir     Path to an output directory (required).
"""
    exit 0
  }


  if (paramExists("input") == paramExists("csv")) {
    exit 1, "ERROR: Please provide either an --input parameter or a --csv parameter"
  }
  
  assertParamExists("reference", "a Cell Ranger reference directory")
  assertParamExists("publishDir", "where output files will be published")

  if (paramExists("csv")) {
    input_ch = Channel.fromPath(params.csv)
      | splitCsv(header: true, sep: ",")
  } else {
    input_ch = Channel.value( params.subMap(["id", "input"]) )
  }

  def reference = file(params.reference)

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

      // process id
      if (li.containsKey("id") && li.id) {
        id_value = li.id
      } else if (!paramExists("csv")) {
        id_value = "run"
      } else {
        exit 1, "ERROR: The provided csv file should contain an 'id' column"
      }
      [ id_value, [ input: input_path, reference: reference] ]
    }
    | view { "Input: $it" }
    | run_wf
    | publish.run(
      map: { [ it[0], [ input: it[1], output: "${it[0]}.h5mu" ] ] },
      auto: [ publish: true ]
    )
    | view { "Output: ${params.publishDir}/${it[1].name}" }
}

/* Cell Ranger Mapping - common workflow
 * 
 * consumed params:
 *   id:                            a sample id for one or more fastq files
 *   data:
 *     input:                       one or more fastq paths, separated with semicolons, paths may be globs
 *     reference:                   Path to a Cell Ranger reference.
 *   output                         a publish dir for the output h5ad files
 * output format:               [ id, h5ad, params ]
 *   value id:                      a sample id for one or more fastq files
 *   value output:                  directory of mapped fastq reads
 *   value params:                  the params object, which may already have sample specific overrides
 */
workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | cellranger_count
    | cellranger_count_split
    | from_10xh5_to_h5mu.run(
      mapData: { it.filtered_h5 }
    )

  emit:
  output_ch
}


/* Cell Ranger Mapping - Integration testing
 */
workflow test_wf {
  
  output_ch =
    Channel.value(
      [
        "foo",
        [
          input: file(params.rootDir + "/resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq"),
          reference: file(params.rootDir + "/resources_test/cellranger_tiny_fastq/cellranger_tiny_ref")
        ],
        params
      ]
    )
    | view { "Input: [${it[0]}, ${it[1]}, params]" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain three elements; [id, file, params]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[1]}"
      "Output: [${output[0]}, ${output[1]}, params]"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    // | check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}