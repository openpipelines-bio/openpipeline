nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { merge } from targetDir + '/dataflow/merge/main.nf'
include { publish }  from targetDir + '/transfer/publish/main.nf'

include { split_modalities_workflow as split_modalities_workflow } from workflowDir + '/multiomics/full_pipeline/main.nf'
include { multisample_processing_workflow as multisample_processing_workflow } from workflowDir + '/multiomics/full_pipeline/main.nf'
include { integration_setup_workflow as integration_setup_workflow } from workflowDir + '/multiomics/full_pipeline/main.nf'

include { readConfig; helpMessage; readCsv; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap } from workflowDir + "/utils/DataflowHelper.nf"
config = readConfig("$workflowDir/multiomics/multisample/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf

}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | preprocessInputs("config": config)
    | setWorkflowArguments (
        "split_modalities_args": [:],
        "rna_multisample_args": [
          "filter_with_hvg_var_output": "filter_with_hvg_var_output",
          "filter_with_hvg_obs_batch_key": "filter_with_hvg_obs_batch_key",
          "var_qc_metrics": "var_qc_metrics",
          "top_n_vars": "top_n_vars"
        ],
        "prot_multisample_args": [:],
          "integration_args_rna": [
            "var_pca_feature_selection": "filter_with_hvg_var_output", // run PCA on highly variable genes only
            "pca_overwrite": "pca_overwrite"
          ],
        "integration_args_prot": ["pca_overwrite": "pca_overwrite"],
        "publish": ["output": "output"],
    )
    | getWorkflowArguments(key: "split_modalities_args")
    | split_modalities_workflow
    | pmap {id, data ->
        def new_data = ["input": data.output]
        [id, new_data]
    }
    | multisample_processing_workflow
    | groupTuple(by: 0, sort: "hash")
    | map { list -> 
        def id = list[0]
        def new_input = list[1].collect({it.output})
        def other_arguments = list[2][0] // Get the first set, the other ones are copies
        def modalities_list = ["modalities": list[-1].collect({it.modality}).unique()]
        [id, ["input": new_input]] + other_arguments + modalities_list
      }
    | merge.run(args: [ output_compression: "gzip" ])
    | pmap {id, data ->
        def new_data = ["input": data.output]
        [id, new_data]
    }
    | integration_setup_workflow
    | pmap {id, data ->
        def new_data = ["input": data.output]
        [id, new_data]
    }
    | getWorkflowArguments(key: "publish")
    | publish
    | map {list -> [list[0], list[1]]}

  emit:
  output_ch
}

workflow test_wf {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  
  testParams = [
      param_list: [
        [
          id: "test",
          input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
          publish_dir: "foo/"
        ],
        [
          id: "test2",
          input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
          publish_dir: "foo/"
        ],
      ]
    ]


  output_ch =
    channelFromParams(testParams, config)
      | view { "Input: $it" }
      | run_wf
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        print "output_list: $output_list"
        assert output_list.size() == 2 : "output channel should contain two events"
        assert output_list.collect({it[0]}).sort() == ["test", "test2"] : "First output ID should be 'test'"
      }
  
}

workflow test_wf2 {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  
  testParams = [
      input: params.resources_test + "/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_mms.h5mu",
      pca_overwrite: true,
      id: "test",
      publish_dir: "foo/",
      output: "test.h5mu"
    ]


  output_ch =
    channelFromParams(testParams, config)
      | view { "Input: $it" }
      | run_wf
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file], was $output"
        assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        print "output_list: $output_list"
        assert output_list.size() == 1 : "output channel should contain two events"
        assert output_list.collect({it[0]}).sort() == ["test"] : "First output ID should be 'test'"
      }
  
}