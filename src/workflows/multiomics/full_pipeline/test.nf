nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { remove_modality }  from targetDir + '/filter/remove_modality/main.nf'
include { full_pipeline } from targetDir + "/workflows/multiomics/full_pipeline/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
    [
      id: "mouse",
      input: resources_test.resolve("concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      publish_dir: "foo/",
      rna_min_counts: 2
    ],
    [
      id: "human",
      input: resources_test.resolve("concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      publish_dir: "foo/",
      rna_min_counts: 2
    ]
  ])
  | map{ state -> [state.id, state] }
  | full_pipeline
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
  
}

workflow test_wf2 {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList([
    [
        id: "pbmc",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        var_name_mitochondrial_genes: 'mitochondrial',
        rna_min_counts: 2,
        prot_min_counts: 3,
        add_id_to_obs: true,
        add_id_make_observation_keys_unique: true,
        add_id_obs_output: "sample_id"
    ],
    [
      id: "pbmc_with_more_args",
      input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
      rna_min_counts: 2,
      rna_max_counts: 1000000,
      rna_min_genes_per_cell: 1,
      rna_max_genes_per_cell: 1000000,
      rna_min_cells_per_gene: 1,
      rna_min_fraction_mito: 0.0,
      rna_max_fraction_mito: 1.0,
      prot_min_counts: 3,
      prot_max_counts: 1000000,
      prot_min_proteins_per_cell: 1,
      prot_max_proteins_per_cell: 1000000,
      prot_min_cells_per_protein: 1,
      var_name_mitochondrial_genes: 'mitochondrial',
      obs_name_mitochondrial_fraction: 'fraction_mitochondrial',
      add_id_to_obs: true,
      add_id_make_observation_keys_unique: true,
      add_id_obs_output: "sample_id"
    ],
  ])
  | map{ state -> [state.id, state] }
  | full_pipeline
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    // The result of this pipeline is always 1 merged sample, regardless of the number of input samples. 
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
}

workflow test_wf3 {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  input_ch = Channel.fromList([
      [
        id: "mouse",
        input: resources_test.resolve("concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
        publish_dir: "foo/",
        rna_min_counts: 2,
        var_qc_metrics: "highly_variable",
        filter_with_hvg_var_output: "highly_variable",
      ],
      [
        id: "human",
        input: resources_test.resolve("concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
        publish_dir: "foo/",
        rna_min_counts: 2,
        var_qc_metrics: "highly_variable",
        filter_with_hvg_var_output: "highly_variable",
      ]
    ])
    | map{ state -> [state.id, state] }

    human_ch = input_ch
      | filter{it[0] == "human"}
      | remove_modality.run(
        fromState: { id, state ->
          [
            "input": state.input,
            "modality": "atac"
          ]
        },
        toState: ["input": "output"]
      )

    mouse_ch = input_ch
      | filter{it[0] == "mouse"}
      | remove_modality.run(
        fromState: { id, state ->
          [
            "input": state.input,
            "modality": "rna"
          ]
        },
        toState: ["input": "output"]
      )

    output_ch_test_2 = human_ch.concat(mouse_ch)
      | full_pipeline
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
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