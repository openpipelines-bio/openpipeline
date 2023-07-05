nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { split_modalities } from targetDir + '/dataflow/split_modalities/main.nf'
include { merge } from targetDir + '/dataflow/merge/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'
include { run_wf as rna_multisample } from workflowDir + '/multiomics/rna_multisample/main.nf'
include { run_wf as prot_multisample } from workflowDir + '/multiomics/prot_multisample/main.nf'
include { run_wf as initialize_integration_rna } from workflowDir + '/multiomics/integration/initialize_integration/main.nf'
include { run_wf as initialize_integration_prot } from workflowDir + '/multiomics/integration/initialize_integration/main.nf'
include { splitStub } from workflowDir + '/multiomics/multisample/split_stub.nf'

include { readConfig; helpMessage; readCsv; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include {  setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap } from workflowDir + "/utils/DataflowHelper.nf"
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

  arguments_parsed_ch = input_ch
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
          "output": "output"
        ],
        "integration_args_prot": ["output": "output"]
    )
    | getWorkflowArguments(key: "split_modalities_args")

  split_ch = arguments_parsed_ch
    | filter{!workflow.stubRun}
    | split_modalities

  split_stub_ch = arguments_parsed_ch
    | filter{workflow.stubRun}
    | map {it -> [it[0], it[1].input, it[2]]}
    | splitStub
    | map {it -> [it[0], ["output": it[1], "output_types": it[2]], it[3]]}

  start_ch = split_ch.concat(split_stub_ch)
    // combine output types csv
    | pFlatMap {id, data, passthrough ->
      def outputDir = data.output
      def types = readCsv(data.output_types.toString())
      
      types.collect{ dat ->
        // def new_id = id + "_" + dat.name
        def new_id = id // it's okay because the channel will get split up anyways
        def new_data = outputDir.resolve(dat.filename)
        def new_passthrough = passthrough
        [ new_id, new_data, new_passthrough, [ modality: dat.name ]]
      }
    }

  modality_processors = [
    ["id": "rna", "processor": rna_multisample],
    ["id": "prot", "processor": prot_multisample]
  ]
  known_modalities = modality_processors.collect{it.id}

  mod_chs = modality_processors.collect{ modality_processor ->    
    // Reformat arguments to serve to the multisample processing
    input_ms_ch = start_ch
      // Select the files corresponding to the currently selected modality
      | filter{ it[3].modality == modality_processor.id }
      | view {"Start ch: $it"}
    // Run the multisample processing if defined, otherwise just concatenate samples together
    out_ch = (
      modality_processor.processor ? \
        input_ms_ch
          | getWorkflowArguments(key: ("$modality_processor.id" + "_multisample_args").toString() )
          | pmap {id, input_args, other_arguments -> 
              [id, ["add_id_to_obs": false, "sample_id": [id]] + input_args, other_arguments]
          } 
          | view { "input multichannel-$modality_processor.id: $it" }
          | modality_processor.processor : \
        input_ms_ch
    )
    return out_ch
  }
    
  // Keep unknown modalities as well
  unknown_channel = start_ch
    | filter { ! known_modalities.contains(it[3].modality.toString())}


  // Concat channel if more than one modality was found
  merge_ch = unknown_channel.concat(*mod_chs)
  for_integration_ch = merge_ch
    | groupTuple(by: 0)
    | map { list -> 
      def id = list[0]
      def new_input = ["input": list[1]]
      def other_arguments = list[2][0] // Get the first set, the other ones are copies
      def modalities = list[3]
      def modalities_list = ["modalities": modalities.collect{it.modality}]
      [id, new_input, other_arguments, modalities_list]
    }
    | merge.run(args: [ output_compression: "gzip" ])
  
  integration_processors = [
    [id: "rna", "workflow": initialize_integration_rna, "args": ["layer": "log_normalized", "modality": "rna"]],
    [id: "prot", "workflow": initialize_integration_prot, "args": ["layer": "clr", "modality": "prot"]],
  ]

  output_ch = integration_processors.inject(for_integration_ch){ channel_in, processor ->
    channel_out_integrated = channel_in
      | filter{it[3].modalities.contains(processor.id)}
      | pmap {id, input_args -> [id, ["input": input_args] + processor.args]}
      | processor.workflow
    ch_in_unmodified = channel_in
      | filter{ !(it[3].modalities.contains(processor.id)) }
    return channel_out_integrated.concat(ch_in_unmodified)
  }
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
        assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        print "output_list: $output_list"
        assert output_list.size() == 2 : "output channel should contain two events"
        assert output_list.collect({it[0]}).sort() == ["test", "test2"] : "First output ID should be 'test'"
      }
  
}