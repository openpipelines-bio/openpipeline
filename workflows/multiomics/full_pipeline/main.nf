nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { add_id } from targetDir + "/metadata/add_id/main.nf"
include { split_modalities } from targetDir + '/dataflow/split_modalities/main.nf'
include { merge } from targetDir + '/dataflow/merge/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'
include { remove_modality }  from targetDir + '/filter/remove_modality/main.nf'
include { publish }  from targetDir + '/transfer/publish/main.nf'
include { run_wf as rna_singlesample } from workflowDir + '/multiomics/rna_singlesample/main.nf'
include { run_wf as rna_multisample } from workflowDir + '/multiomics/rna_multisample/main.nf'
include { run_wf as prot_singlesample } from workflowDir + '/multiomics/prot_singlesample/main.nf'
include { run_wf as prot_multisample } from workflowDir + '/multiomics/prot_multisample/main.nf'
include { run_wf as initialize_integration_rna } from workflowDir + '/multiomics/integration/initialize_integration/main.nf'
include { run_wf as initialize_integration_prot } from workflowDir + '/multiomics/integration/initialize_integration/main.nf'
include { splitStub } from workflowDir + '/multiomics/full_pipeline/split_stub.nf'

include { readConfig; helpMessage; readCsv; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include {  setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap; strictMap as smap } from workflowDir + "/utils/DataflowHelper.nf"
config = readConfig("$workflowDir/multiomics/full_pipeline/config.vsh.yaml")

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
      | add_arguments
      | getWorkflowArguments(key: "add_id_args")
      | add_id_workflow
      | getWorkflowArguments(key: "split_modalities_args")
      | split_modalities_workflow
      | singlesample_processing_workflow
      | concat_workflow
      | multisample_processing_workflow
      | toSortedList{ a, b -> b[0] <=> a[0] }
      | map { list -> 
          def new_input = list.collect{it[1].input}
          def other_arguments = list[0][2] // Get the first set, the other ones are copies
          def modalities_list = ["modalities": list.collect({it[-1].modality}).unique()]
          ["merged", ["input": new_input]] + other_arguments + modalities_list
      }
      | merge.run(args: [ output_compression: "gzip" ])
      | integration_setup_workflow
      | getWorkflowArguments(key: "publish")
      | publish
      | map {list -> [list[0], list[1]]}
      

  emit:
    output_ch
}

// =================================
// === start of helper workflows ===
// =================================

workflow add_arguments {
  // Process the specified arguments for the workflow run and perform basic assertions on argument sanity.
  take:
    input_ch

  main:
    parsed_arguments_ch = input_ch
      | preprocessInputs("config": config)
      | pmap {id, args ->
        def var_qc_default = [args.filter_with_hvg_var_output]
        if (args.var_name_mitochondrial_genes) {
          var_qc_default.add(args.var_name_mitochondrial_genes)
        }
        def new_args = args + ["var_qc_metrics": var_qc_default.join(",")]
        [id, new_args]
      }
      | setWorkflowArguments (
          "add_id_args": ["input": "input",
                          "make_observation_keys_unique": "add_id_make_observation_keys_unique",
                          "obs_output": "add_id_obs_output",
                          "add_id_to_obs": "add_id_to_obs"],
          "split_modalities_args": [:],
          "rna_singlesample_args": [
            "min_counts": "rna_min_counts",
            "max_counts": "rna_max_counts",
            "min_genes_per_cell": "rna_min_genes_per_cell",
            "max_genes_per_cell": "rna_max_genes_per_cell",
            "min_cells_per_gene": "rna_min_cells_per_gene",
            "min_fraction_mito": "rna_min_fraction_mito",
            "max_fraction_mito": "rna_max_fraction_mito",
            "var_name_mitochondrial_genes": "var_name_mitochondrial_genes",
            "var_gene_names": "var_gene_names",
            "mitochondrial_gene_regex": "mitochondrial_gene_regex"
          ],
          "prot_singlesample_args": [
            "min_counts": "prot_min_counts",
            "max_counts": "prot_max_counts",
            "min_proteins_per_cell": "prot_min_proteins_per_cell",
            "max_proteins_per_cell": "prot_max_proteins_per_cell",
            "min_cells_per_protein": "prot_min_cells_per_protein"
          ],
          "rna_multisample_args": [
            "filter_with_hvg_var_output": "filter_with_hvg_var_output",
            "filter_with_hvg_obs_batch_key": "filter_with_hvg_obs_batch_key",
            "var_qc_metrics": "var_qc_metrics",
            "top_n_vars": "top_n_vars"
          ],
          "prot_multisample_args": [:],
          "integration_args_rna": [
            "var_pca_feature_selection": "filter_with_hvg_var_output", // run PCA on highly variable genes only
          ],
          "integration_args_prot": [:],
          "publish": ["output": "output"]
      )
    
    parsed_arguments_ch 
      | toSortedList
      | map { list ->
        found_output_files = list.collect{it[2].get('publish').getOrDefault("output", null)}.unique()
        assert found_output_files.size() < 2, "The specified output file is not the same for all samples. Found: $found_output_files"
      }


  emit:
    parsed_arguments_ch
}

workflow add_id_workflow {
  // If requested, add the id of the events (samples) a column in .obs. 
  // Also allows to make .obs_names (the .obs index) unique, by prefixing the values with an unique id per .h5mu file.
  // The latter is useful to avoid duplicate observations during concatenation.
  take: 
    parsed_arguments_ch

  main:
    id_added_ch = parsed_arguments_ch
      | filter{ it[1].add_id_to_obs }
      // add ids 
      | pmap { id, data ->
          def new_data = data + [input_id: id]
          [id, new_data]
      }
      | add_id

    no_id_added_ch = parsed_arguments_ch
      | filter{ ! (it[1].add_id_to_obs) }

    samples_with_id_ch = id_added_ch.mix(no_id_added_ch)

  emit:
    samples_with_id_ch
  }

workflow split_modalities_workflow {
  // Split multimodal MuData files into several unimodal MuData files.
  take:
    samples_with_id_ch

  main:
    split_ch = samples_with_id_ch
      | filter{!workflow.stubRun}
      | split_modalities

    split_stub_ch = samples_with_id_ch
      | filter{workflow.stubRun}
      | map {it -> [it[0], it[1].input, it[2]]}
      | splitStub
      | map {it -> [it[0], ["output": it[1], "output_types": it[2]], it[3]]}

    start_ch = split_ch.concat(split_stub_ch)
      // combine output types csv
      | pFlatMap {id, data, passthrough ->
        def outputDir = data.output
        def types = readCsv(data.output_types.toUriString())
        
        types.collect{ dat ->
          // def new_id = id + "_" + dat.name
          def new_id = id // it's okay because the channel will get split up anyways
          def new_data = outputDir.resolve(dat.filename)
          def new_passthrough = passthrough
          [ new_id, ["input": new_data], new_passthrough, [ modality: dat.name ]]
        }
      }

  emit:
    start_ch
}

workflow singlesample_processing_workflow {
  // Perform processing of unimodal single-sample MuData object
  // Each modality and each sample is processed individually.
  take:
    start_ch

  main:
    modality_processors = [
      ["id": "rna", "singlesample": rna_singlesample],
      ["id": "prot", "singlesample": prot_singlesample]
    ]

    known_modalities = modality_processors.collect{it.id}

    mod_chs = modality_processors.collect{ modality_processor ->
      start_ch
        // Select the files corresponding to the currently selected modality
        | filter{ it[3].modality == modality_processor.id }
        | getWorkflowArguments(key: ("$modality_processor.id" + "_singlesample_args").toString() )
        | view { "single-sample-input-$modality_processor.id: $it" }
        | modality_processor.singlesample
    }

    unknown_mods = start_ch
      | filter { ! known_modalities.contains(it[3].modality.toString())}
      | smap { id, args, state, modality ->
        [id, args.input, state, modality]
      }
  
    singlesample_ch = unknown_mods.concat(*mod_chs)

  emit:
    singlesample_ch
}

workflow concat_workflow {
  // Concatenate multiple single-sample unimodal MuData objects back into several multi-sample files.
  // One multi-sample MuData file is created per modality.
  take:
    singlesample_ch

  main:
    concat_ch = singlesample_ch
      | map { lst -> // Put modality name in first element so that we can group on it
        [lst[3].modality] + lst
      }
      | groupTuple(by: 0, sort: "hash") 
      | map { grouped_lst ->
        def new_id = "combined_${grouped_lst[0]}"
        def new_data = [
          "input": grouped_lst[2],
          "input_id": grouped_lst[1]
        ]
        // passthrough is copied, just pick the first
        def new_passthrough = grouped_lst[3][0]
        def new_modalities = grouped_lst[4][0]
        [new_id, new_data] + new_passthrough + new_modalities
      }
      | concat.run(
          // The Ids have already been added in this pipeline
          args: [ add_id_to_obs: false ],
          auto: [ simplifyOutput: false ]
      )
      | smap {id, args, state, modality -> [id, ["input": args.output], state, modality]}
      | view{"Concat output: $it"}

  emit:
    concat_ch

}

workflow multisample_processing_workflow {
  // Perform processing of multisample unimodal MuData files.
  // Each modality is processed individually.
  take:
    concat_ch

  main:
    modality_processors = [
      ["id": "rna", "multisample": rna_multisample],
      ["id": "prot", "multisample": prot_multisample]
    ]
    known_modalities = modality_processors.collect{it.id}
    
    mod_chs = modality_processors.collect{ modality_processor ->
      concat_ch
        | filter{ it[3].modality == modality_processor.id }
        | getWorkflowArguments(key: ("$modality_processor.id" + "_multisample_args").toString() )
        | view { "input multichannel-$modality_processor.id: $it" }
        | modality_processor.multisample
        | smap {id, args, passthrough, modality -> [id, ["input": args], passthrough, modality]}

        | view { "output multichannel-$modality_processor.id: $it" }
    }

    unknown_channel = concat_ch
      | filter { ! known_modalities.contains(it[3].modality.toString())}
      | view { "output multichannel-" + it[3].modality + ": $it" }

    multisample_ch = unknown_channel.concat(*mod_chs)

  emit:
    multisample_ch
}

workflow integration_setup_workflow {
  // Processing of multi-modal multisample MuData files.
  // Performs calculations on samples that have *not* been integrated,
  // and can be considered a "no-integration" workflow.
  take:
    for_integration_ch

  main:
    integration_processors = [
      [id: "rna", "workflow": initialize_integration_rna, "extra_args": ["layer": "log_normalized", "modality": "rna"]],
      [id: "prot", "workflow": initialize_integration_prot, "extra_args": ["layer": "clr", "modality": "prot"]],
    ]

    output_ch = integration_processors.inject(for_integration_ch){ channel_in, processor ->
      channel_out_integrated = channel_in
        | filter{it[3].modalities.contains(processor.id)}
        | getWorkflowArguments(key: ("integration_args_" + processor.id).toString() )
        | pmap {id, input_args -> [id, input_args + processor.extra_args]}
        | view {"integration-input-$processor.id: $it"}
        | processor.workflow
      ch_in_unmodified = channel_in
        | filter{ !(it[3].modalities.contains(processor.id)) }
      return channel_out_integrated.concat(ch_in_unmodified)
    }

  emit:
    output_ch
}

// ===============================
// === start of test workflows ===
// ===============================

workflow test_wf {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
        [
          id: "mouse",
          input: params.resources_test + "/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
          publish_dir: "foo/"
        ],
        [
          id: "human",
          input: params.resources_test + "/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
          publish_dir: "foo/"
        ]
      ],
      rna_min_counts: 2,
      output: "full_pipeline_output.h5mu"
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
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
        assert output_list[0][1].getFileName().toString() == "full_pipeline_output.h5mu" : "Output file should be named 'full_pipeline_output.h5mu'. Found: ${output_list[0][1]}"
      }
  
}


workflow test_wf3 {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
      [
        id: "mouse",
        input: params.resources_test + "/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
        publish_dir: "foo/"
      ],
      [
        id: "human",
        input: params.resources_test + "/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
        publish_dir: "foo/"
      ]
    ],
    rna_min_counts: 2,
    var_qc_metrics: "highly_variable",
    filter_with_hvg_var_output: "highly_variable",
  ]

  input_ch = channelFromParams(testParams, config)
    // store output value in 3rd slot for later use
    // and transform for concat component
    | map { tup ->
      data = tup[1]
      new_data = [ id: data.id, input: data.input ]
      [tup[0], new_data, data] + tup.drop(2)
    }

    human_ch = input_ch
      | filter{it[0] == "human"}
      | remove_modality.run(
        args: [ modality: "atac" ]
      )

    mouse_ch = input_ch
      | filter{it[0] == "mouse"}
      | remove_modality.run(
        args: [ modality: "rna" ]
      )

    output_ch_test_2 = human_ch.concat(mouse_ch)
      // Put back values for other parameters after removing modalities
      | map { tup ->
        tup[2].remove("input")
        new_data = [input: tup[1]]
        new_data.putAll(tup[2])
        [tup[0], new_data]
      }
      | view { "Input test 2: $it" }
      | run_wf
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
      }
}

workflow test_wf2 {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  testParams = [
    param_list: [
        [
          id: "pbmc",
          input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
          var_name_mitochondrial_genes: 'mitochondrial'
        ],
        [
          id: "pbmc_with_more_args",
          input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
          rna_min_counts: 1,
          rna_max_counts: 1000000,
          rna_min_genes_per_cell: 1,
          rna_max_genes_per_cell: 1000000,
          rna_min_cells_per_gene: 1,
          rna_min_fraction_mito: 0,
          rna_max_fraction_mito: 1,
          prot_min_counts: 1,
          prot_max_counts: 1000000,
          prot_min_proteins_per_cell: 1,
          prot_max_proteins_per_cell: 1000000,
          prot_min_cells_per_protein: 1,
          var_name_mitochondrial_genes: 'mitochondrial',
        ],
      ],
      rna_min_counts: 2,
      prot_min_counts: 3
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
        // The result of this pipeline is always 1 merged sample, regardless of the number of input samples. 
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
      }
}

// The following test is supposed to fail. It is used to test the error handling of the pipeline.
// However, there is not way to catch the error originating from the workflow

// workflow test_wf4 {

//   helpMessage(config)

//   // allow changing the resources_test dir
//   params.resources_test = params.rootDir + "/resources_test"
  

//   testParams = [
//     param_list: [
//         [
//           id: "pbmc",
//           input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
//           output: "foo.h5mu"
//         ],
//         [
//           id: "pbmc_with_more_args",
//           input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
//           output: "foo2.h5mu"
//         ],
//       ],

//     ]

//     output_ch = channelFromParams(testParams, config)
//       | run_wf
// }