nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { cellbender_remove_background } from targetDir + "/correction/cellbender_remove_background/main.nf"
include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"
include { from_10xh5_to_h5mu } from targetDir + "/convert/from_10xh5_to_h5mu/main.nf"
include { subset_h5mu } from targetDir + "/filter/subset_h5mu/main.nf"
include { publish } from targetDir + "/transfer/publish/main.nf"

include { readConfig; channelFromParams; preprocessInputs; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/ingestion/cellranger_postprocessing/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  main:

  mid0 = input_ch
    | preprocessInputs("config": config)
    // split params for downstream components
    | setWorkflowArguments(
      correction: [
        "perform_correction": "perform_correction"
      ],
      filter_with_counts: [
        "min_genes": "min_genes",
        "min_counts": "min_counts",
      ],
      publish: [
        "output": "output"
      ]
    )


  // perform correction if so desired
  mid1_corrected = mid0
    | filter{ it[2].correction.perform_correction }
    | cellbender_remove_background
    | pmap{ id, file -> [ id, [ input: file, layer: "corrected" ]]}
    // todo: allow setting the layer
  mid1_uncorrected = mid0
    | filter{ ! it[2].correction.perform_correction }
  mid1 = mid1_corrected.mix(mid1_uncorrected)

  // perform filtering if so desired
  // todo: set layer
  mid2_filtered = mid1
    | filter{ it[2].filter_with_counts.min_genes != null || it[2].filter_with_counts.min_counts != null }
    | getWorkflowArguments(key: "filter_with_counts")
    | filter_with_counts.run(args: [do_subset: true])
  mid2_unfiltered = mid1
    | filter{ it[2].filter_with_counts.min_genes == null && it[2].filter_with_counts.min_counts == null }
  mid2 = mid2_filtered.mix(mid2_unfiltered)
    
  // return output map
  output_ch = mid2
    | pmap { id, data, split_args ->
      file = data instanceof Map ? data.input : data
      [ id, [input: file, output: split_args.publish.output] ]
    }
    | publish.run(auto: [ publish: true, simplifyOutput: false ])

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
      [
        id: "foo",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
        perform_correction: true,
        min_genes: 100,
        min_counts: 1000
      ],
      [
        id: "bar",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
        perform_correction: true
      ],
      [
        id: "zing",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
        min_genes: 100,
        min_counts: 1000
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)

    // first filter and convert to h5mu
    | pmap { id, data -> [ id, [ input: data.input ], data ] }
    | from_10xh5_to_h5mu
    | pmap { id, file -> [ id, [ input: file ] ] }
    | subset_h5mu.run(args: [ number_of_observations: 100000 ])
    | pmap { id, file, orig_params -> [id, orig_params + [ input: file ] ] }

    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 3 : "output channel should contain three events"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}