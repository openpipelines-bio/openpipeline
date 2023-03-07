nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { add_id } from targetDir + "/metadata/add_id/main.nf"
include { split_modalities } from targetDir + '/dataflow/split_modalities/main.nf'
include { merge } from targetDir + '/dataflow/merge/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'
include { remove_modality }  from targetDir + '/filter/remove_modality/main.nf'
include { run_wf as rna_singlesample } from workflowDir + '/multiomics/rna_singlesample/main.nf'
include { run_wf as rna_multisample } from workflowDir + '/multiomics/rna_multisample/main.nf'
include { run_wf as prot_singlesample } from workflowDir + '/multiomics/prot_singlesample/main.nf'
include { run_wf as integration } from workflowDir + '/multiomics/integration/main.nf'

include { readConfig; helpMessage; readCsv; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include {  setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap } from workflowDir + "/utils/DataflowHelper.nf"
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
  start_ch = input_ch
    | preprocessInputs("config": config)
    | setWorkflowArguments (
        "add_id_args": ["input": "input"],
        "split_modalities_args": [:],
        "rna_singlesample_args": [
          "min_counts": "rna_min_counts",
          "max_counts": "rna_max_counts",
          "min_genes_per_cell": "rna_min_vars_per_cell",
          "max_genes_per_cell": "rna_max_vars_per_cell",
          "min_cells_per_gene": "rna_min_cells_per_gene",
          "min_fraction_mito": "rna_min_fraction_mito",
          "max_fraction_mito": "rna_max_fraction_mito",
        ],
        "prot_singlesample_args": [
          "min_counts": "prot_min_counts",
          "max_counts": "prot_max_counts",
          "min_genes_per_cell": "prot_min_vars_per_cell",
          "max_genes_per_cell": "prot_max_vars_per_cell",
          "min_cells_per_gene": "prot_min_cells_per_protein",
          "min_fraction_mito": "prot_min_fraction_mito",
          "max_fraction_mito": "prot_max_fraction_mito",
        ],
        "rna_multisample_args": [
          "filter_with_hvg_var_output": "filter_with_hvg_var_output",
          "filter_with_hvg_obs_batch_key": "filter_with_hvg_obs_batch_key",
          "var_qc_metrics": "var_qc_metrics",
          "top_n_vars": "top_n_vars"
        ],
        "prot_multisample_args": [:],
        "integration_args": [
          "obs_covariates": "obs_covariates",
          "var_pca_feature_selection": "filter_with_hvg_var_output", // run PCA on highly variable genes only
          "rna_theta": "rna_harmony_theta",
          "leiden_resolution": "leiden_resolution",
        ]
    )
    // add ids to obs_names and to .obs[sample_id]
    | getWorkflowArguments(key: "add_id_args")

    | pmap { id, data ->
        def new_data = data + [
          input_id: id,
          make_observation_keys_unique: true,
          obs_output: "sample_id"
        ]
        [id, new_data]
    }
    | add_id

    // split by modality
    | getWorkflowArguments(key: "split_modalities_args")
    | split_modalities

    // combine output types csv
    | pFlatMap {id, data, passthrough ->
      def outputDir = data.output
      def types = readCsv(data.output_types.toString())
      
      types.collect{ dat ->
        // def new_id = id + "_" + dat.name
        def new_id = id // it's okay because the channel will get split up anyways
        def new_data = outputDir.resolve(dat.filename)
        def new_passthrough = passthrough
        [ new_id, ["input": new_data], new_passthrough, [ modality: dat.name ]]
      }
    }

    modality_processors = [
      ["id": "rna", "singlesample": rna_singlesample, "multisample": rna_multisample],
      ["id": "prot", "singlesample": prot_singlesample, "multisample": null]
    ]
    known_modalities = modality_processors.collect{it.id}

    mod_chs = modality_processors.collect{ modality_processor ->
      // Select the files corresponding to the currently selected modality
      mod_ch = start_ch
        | filter{ it[3].modality == modality_processor.id }
        | getWorkflowArguments(key: ("$modality_processor.id" + "_singlesample_args").toString() )
        | view { "single-sample-input-$modality_processor.id: $it" }
      // Run the single-sample processing if defined
      ss_ch = (modality_processor.singlesample ? \
               mod_ch | modality_processor.singlesample : \
               mod_ch)
      
      // Reformat arguments to serve to the multisample processing
      input_ms_ch = ss_ch
        | view { "single-sample-output-$modality_processor.id: $it" }
        | toSortedList{ a, b -> b[0] <=> a[0] }
        | filter { it.size() != 0 } // filter when event is empty
        | map{ list -> 
          def new_data = ["sample_id": list.collect{it[0]}, "input": list.collect{it[1]}]
          ["combined_$modality_processor.id", new_data] + list[0].drop(2)
        }
        | getWorkflowArguments(key: ("$modality_processor.id" + "_multisample_args").toString() )
        | view { "input multichannel-$modality_processor.id: $it" }
      
      // Run the multisample processing if defined, otherwise just concatenate samples together
      out_ch = (
        modality_processor.multisample ? \
          input_ms_ch | modality_processor.multisample : \
          input_ms_ch | concat.run(
            key: "concat_" + modality_processor.id,
            renameKeys: [input_id: "sample_id"]
          )
      )
      return out_ch
    }
    
  // Keep and concat unknown modalities as well
  unknown_channel = start_ch
    | filter { ! known_modalities.contains(it[3].modality.toString())}
    | map { lst -> // Put modality name in first element so that we can group on it
        [lst[3].modality] + lst
      }
    | groupTuple(by: 0)
    // [ String, List[String], List[Map[String, Any]], ... ]
    | map { grouped_lst ->
      def new_id = "combined_${grouped_lst[0]}"
      def new_data = [
        "input": grouped_lst[2].collect{it.input},
        "input_id": grouped_lst[1]
      ]
      // passthrough is copied, just pick the first
      def new_passthrough = grouped_lst.drop(3).collect{it[0]}
      [new_id, new_data] + new_passthrough
    }
    | concat

  // Concat channel if more than one modality was found
  merge_ch = unknown_channel.concat(*mod_chs)
  output_ch = merge_ch
    | toSortedList{ a, b -> b[0] <=> a[0] }
    | map { list -> 
      ["merged", list.collect{it[1]}] + list[0].drop(2)
    }
    | merge
    | getWorkflowArguments(key: "integration_args")
    | pmap {id, input_args -> [id, input_args + [layer: "log_normalized"]]}
    | integration
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
      obs_covariates: "sample_id",
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
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
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
    obs_covariates: [],
    rna_min_counts: 2,
    var_qc_metrics: "highly_variable",
    filter_with_hvg_var_output: "highly_variable",
    rna_harmony_theta: 3,
    leiden_resolution: 2,
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
          publish_dir: "test2/"
        ],
      ],
      obs_covariates: "sample_id",
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
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
      }
}