nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf" params(params)
include { filter_with_scrublet } from targetDir + "/filter/filter_with_scrublet/main.nf" params(params)
include { do_filter } from targetDir + "/filter/do_filter/main.nf" params(params)
include { do_filter as do_filter2 } from targetDir + "/filter/do_filter/main.nf" params(params)
include { lognorm } from targetDir + '/normalize/lognorm/main.nf' params(params)
include { filter_with_hvg } from targetDir + '/filter/filter_with_hvg/main.nf' params(params)
include { pca } from targetDir + '/dimred/pca/main.nf' params(params)
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf' params(params)
include { umap } from targetDir + '/dimred/umap/main.nf' params(params)
include { leiden } from targetDir + '/cluster/leiden/main.nf' params(params)

include { publish } from targetDir + "/transfer/publish/main.nf" params(params)
include { overrideOptionValue; has_param; check_required_param } from workflowDir + "/utils/utils.nf" params(params)


workflow {
  if (has_param("help")) {
    log.info """TX Processing - CLI workflow

A workflow for running the default RNA processing components.
This workflow can be run on a single input or in batch, see below.

Parameters (Single input mode):
  --id       ID of the sample (optional).
  --input    Path to the sample (required).
  --output   Path to an output directory (required).
  
Parameters (Batch mode):
  --csv      A csv file containing columns 'id' and 'input' (required).
  --output   Path to an output directory (required).
"""
    exit 0
  }


  if (has_param("input") == has_param("csv")) {
    exit 1, "ERROR: Please provide either an --input parameter or a --csv parameter"
  }
  
  check_required_param("output", "where output files will be published")

  def multirun = has_param("csv")
  if (multirun) {
    input_ch = Channel.fromPath(params.csv)
      | splitCsv(header: true, sep: ",")
  } else {
    input_ch = Channel.value( params.subMap(["id", "input"]) )
  }

  input_ch
    | map { li ->
      // process input
      if (li.containsKey("input") && li.input) {
        input_path = li.input.split(";").collect { path -> file(path) }.flatten()
      } else {
        exit 1, multirun ? 
          "ERROR: The provided csv file should contain an 'input' column" : 
          "ERROR: Please specify an '--input' parameter"
      }

      // process id
      if (li.containsKey("id") && li.id) {
        id_value = li.id
      } else if (!multirun) {
        id_value = "run"
      } else {
        exit 1, "ERROR: The provided csv file should contain an 'id' column"
      }
      [ id_value, [ input: input_path ], params ]
    }
    | view { "before run_wf: ${it[0]} - ${it[1]}" }
    | run_wf
    | view { "after run_wf: ${it[0]} - ${it[1]}" }
    | map { overrideOptionValue(it, "publish", "output", "${params.output}/${it[0]}.h5mu") }
    | publish
}

/*
TX Processing - Common workflow

A workflow for running the default RNA processing components.

input channel event format: [ id, file, params ]
  value id:                      an event id
  value file:                    an h5mu input file
  value params:                  the params object, which may already have sample specific overrides
output channel event format: [ id, file, params ]
  value id:                      same as input
  value file:                    an h5mu output file
  value params:                  same as input params
*/
workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | filter_with_counts
    | filter_with_scrublet
    | map { overrideOptionValue(it, "do_filter", "obs_filter", "filter_with_counts:filter_with_scrublet") }
    | map { overrideOptionValue(it, "do_filter", "var_filter", "filter_with_counts") }
    | do_filter
    | lognorm
    | filter_with_hvg
    | map { overrideOptionValue(it, "do_filter", "obs_filter", "") }
    | map { overrideOptionValue(it, "do_filter", "var_filter", "filter_with_hvg") }
    | do_filter2
    | pca
    | find_neighbors
    | leiden
    | umap

  emit:
  output_ch
}


/*
TX Processing - Integration testing

A workflow for running the default RNA processing components.
*/
workflow test_wf {
  
  output_ch =
    Channel.value(
      [
        "foo",
        file(params.rootDir + "/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        params
      ]
    )
    | view { "Input: [${it[0]}, ${it[1]}, params]" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain three elements; [id, file, params]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[0][1]}"
      "Output: [${output[0]}, ${output[1]}, params]"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}