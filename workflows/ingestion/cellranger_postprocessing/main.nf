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

  // perform correction if so desired
  mid1_corrected = mid0
    | filter{ it[1].perform_correction }
    | cellbender_remove_background.run(
      fromState: { id, state ->
        [
          input: state.input,
          epochs: state.cellbender_epochs,
          output_layer: "cellbender_corrected",
          output_compression: "gzip"
        ]
      },
      toState: { id, output, state -> 
        state + [input: output.output, layer: "cellbender_corrected"]
      }
    )
  mid1_uncorrected = mid0
    | filter{ ! it[1].perform_correction }
  mid1 = mid1_corrected.mix(mid1_uncorrected)

  // perform filtering if so desired
  mid2_filtered = mid1
    | filter{ it[1].min_genes != null || it[1].min_counts != null }
    | filter_with_counts.run(
      fromState: { id, state ->
        [
          input: state.input,
          min_genes: state.min_genes,
          min_counts: state.min_counts,
          layer: state.layer,
          output_compression: "gzip",
          do_subset: true
        ]
      },
      toState: [input: "output"]
    )
  mid2_unfiltered = mid1
    | filter{ it[1].min_genes == null && it[1].min_counts == null }
  mid2 = mid2_filtered.mix(mid2_unfiltered)
    
  // return output map
  output_ch = mid2
    | publish.run(
      fromState: [ input: "input", output: "output" ],
      auto: [ publish: true ]
    )

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
        min_counts: 1000,
        cellbender_epochs: 5
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)

    // first filter and convert to h5mu
    | from_10xh5_to_h5mu.run(
      fromState: ["input"],
      toState: ["input": "output"]
    )

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
      assert output_list.size() == 1 : "output channel should contain one event"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}

workflow test_wf2 {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
      [
        id: "zing",
        input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
        perform_correction: false,
        min_genes: 100,
        min_counts: 1000,
        cellbender_epochs: 5
      ]
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    
    // first filter and convert to h5mu
    | from_10xh5_to_h5mu.run(
      fromState: ["input"],
      toState: ["input": "output"]
    )

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
      assert output_list.size() == 1 : "output channel should contain one event"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}