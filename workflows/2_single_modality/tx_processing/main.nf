nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf" params(params)
include { filter_with_scrublet } from targetDir + "/filter/filter_with_scrublet/main.nf" params(params)
include { lognorm } from targetDir + '/normalize/lognorm/main.nf' params(params)
include { hvg_scanpy } from targetDir + '/hvg/hvg_scanpy/main.nf' params(params)
include { pca } from targetDir + '/dimred/pca/main.nf' params(params)
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf' params(params)
include { umap } from targetDir + '/dimred/umap/main.nf' params(params)
include { leiden } from targetDir + '/cluster/leiden/main.nf' params(params)

include { publish } from targetDir + "/transfer/publish/main.nf" params(params)
include { overrideOptionValue; has_param; check_required_param } from workflowDir + "/utils/utils.nf" params(params)

/*
TX Processing - CLI workflow

A workflow for running the default RNA processing components.
Exactly one of '--input' and '--csv' must be passed as a parameter.

Parameters:
  --id       ID of the sample, optional.
  --input    Path to the sample.
  --csv      A CSV file with required columns 'input' and optional columns 'id'.
  --output   Path to an output directory.
*/
workflow {
  if (has_param("help")) {
    log.info """TX Processing - CLI workflow

A workflow for running the default RNA processing components.
Exactly one of '--input' and '--csv' must be passed as a parameter.

Parameters:
  --id       ID of the sample, optional.
  --input    Path to the sample.
  --csv      A CSV file with required columns 'input' and optional columns 'id'.
  --output   Path to an output directory."""
    exit 0
  }

  if (has_param("input") == has_param("csv")) {
    exit 1, "ERROR: Please provide either an --input parameter or a --csv parameter"
  }
  if (has_param("input")) {
    if (has_param("id")) {
      input_ch = Channel.value( [ params.id, file(params.input), params ])
    } else {
      input_ch = 
        Channel.fromPath(params.input)
        | map { input_file -> [ input_file.baseName, input_file, params ]}
    }
  } else if (has_param("csv")) {
    input_ch = 
      Channel.fromPath(params.csv)
      | splitCsv(header: true, sep: ",")
      | map { li -> 
        if (!li.containsKey("input")) {
          exit 1, "ERROR: The provided csv file should contain an 'input' column"
        }
        input_path = file(li.input)
        // todo: check if input_path has length 1
        if (li.containsKey("id")) {
          id = li.id
        } else {
          // derive pathname from input
          id = input_path.baseName
        }
        [ id, input_path, params ] 
      }
  }
  
  check_required_param("output", "where output files will be published")

  input_ch
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
    | lognorm
    | hvg_scanpy
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
        params.rootDir + "/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
        params
      ]
    )
    | run_wf
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
  
  output_list = output_ch.toList()
  /* output_list = [ 
    [id, file, params]
  ]*/

  print(output_list)
  System.out.flush()

  assert output_list.size() == 1 : "output_ch should contain one element"
  assert output_list.every{it.size() == 3} : "events in output_ch should contain three elements; [id, file, params]"
  assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
  assert output_list[0][1].matches('.*\\.h5mu\$') : "Output file should be a h5mu file"
}