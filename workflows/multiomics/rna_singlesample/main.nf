nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"
include { filter_with_scrublet } from targetDir + "/filter/filter_with_scrublet/main.nf"
include { do_filter } from targetDir + "/filter/do_filter/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/multiomics/rna_singlesample/config.vsh.yaml")

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }

}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // store output value in 3rd slot for later use
    | map { tup -> 
      data = tup[1]
      if (data !instanceof Map) {
        data = [ input: data ]
      }
      [ tup[0], data, data.clone() ] + tup.drop(2)
    }

    // cell filtering
    | filter_with_counts
    | do_filter.run(
      args: [ obs_filter: "filter_with_counts" ]
    )

    // retrieve output value
    | map { tup -> 
      orig_data = tup[2]
      [ tup[0], [ input: tup[1] ] + (orig_data instanceof Map ? orig_data.subMap("output") : [ output: orig_data ] )] + tup.drop(3)
    }

    // doublet calling
    | filter_with_scrublet.run(
      auto: [ publish: true ]
    )
    // TODO: ambient rna correction

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
  ]

  output_ch =
    viashChannel(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}