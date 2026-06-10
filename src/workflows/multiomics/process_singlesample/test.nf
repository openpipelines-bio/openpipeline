nextflow.enable.dsl=2
targetDir = params.rootDir + "/target"

include { process_singlesample } from targetDir + "/nextflow/workflows/multiomics/process_singlesample/main.nf"
include { workflow_test } from targetDir + "/_test/nextflow/test_workflows/multiomics/process_singlesample/workflow_test/main.nf"
include { workflow2_test } from targetDir + "/_test/nextflow/test_workflows/multiomics/process_singlesample/workflow2_test/main.nf"
include { workflow3_test } from targetDir + "/_test/nextflow/test_workflows/multiomics/process_singlesample/workflow3_test/main.nf"
include { workflow4_test } from targetDir + "/_test/nextflow/test_workflows/multiomics/process_singlesample/workflow4_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  input_ch = Channel.fromList([
    [
      id: "mouse",
      input: resources_test.resolve("concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      publish_dir: "foo/",
      output: "test.h5mu",
    ],
    [
      id: "human",
      input: resources_test.resolve("concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      publish_dir: "foo/",
      output: "test.h5mu",

    ]
  ])
  | map{ state -> [state.id, state] }
  | process_singlesample.run(
    toState: { id, output, state -> output + [orig_input: state.input] }
  )

  assert_ch = input_ch
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith("test.h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 2 : "output channel should contain two events, got ${output_list.size()}"
    assert output_list.collect({it[0]}).sort() == ["human", "mouse"] : "Output id's should contain `mouse` and `human`.'"
  }
      
  test_ch = input_ch
    | workflow_test.run(
      fromState: [
        "input": "output",
        "orig_input": "orig_input"
      ],
    )
}

workflow test_wf2 {

  resources_test = file(params.resources_test)

  input_ch = Channel.fromList([
    [
      id: "pbmc",
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
      var_name_ribosomal_genes: 'ribosomal',
      obs_name_ribosomal_fraction: 'fraction_ribosomal',
      add_id_to_obs: true,
      add_id_make_observation_keys_unique: true,
      do_subset: true,
      add_id_obs_output: "sample_id",
      intersect_obs: true,
      skip_qc_metrics: true,
      output: "pbmc_test.h5mu"
    ]
  ])
  | map{ state -> [state.id, state] }
  | process_singlesample.run(
    toState: { id, output, state -> output + [orig_input: state.input] }
  )

  assert_ch = input_ch
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith("pbmc_test.h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 1 : "output channel should contain one event, got ${output_list.size()}"
    assert output_list[0][0] == "pbmc" : "Output id should be `pbmc`."
  }
      
  test_ch = input_ch
    | workflow2_test.run(
      fromState: [
        "input": "output",
        "orig_input": "orig_input"
      ],
    )
}

workflow test_wf3 {

  resources_test = file(params.resources_test)

  input_ch = Channel.fromList([
    [
      id: "pbmc_scrublet_threshold",
      input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
      rna_min_counts: 2,
      rna_max_counts: 1000000,
      rna_min_genes_per_cell: 1,
      rna_max_genes_per_cell: 1000000,
      rna_min_cells_per_gene: 1,
      prot_min_counts: 3,
      prot_max_counts: 1000000,
      prot_min_proteins_per_cell: 1,
      prot_max_proteins_per_cell: 1000000,
      prot_min_cells_per_protein: 1,
      add_id_to_obs: true,
      add_id_make_observation_keys_unique: true,
      add_id_obs_output: "sample_id",
      scrublet_score_threshold: 0.1,
      scrublet_do_subset: true,
      skip_qc_metrics: true,
      output: "pbmc_scrublet_threshold_test.h5mu"
    ]
  ])
  | map{ state -> [state.id, state] }
  | process_singlesample.run(
    toState: { id, output, state -> output + [orig_input: state.input, scrublet_score_threshold: state.scrublet_score_threshold] }
  )

  assert_ch = input_ch
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith("pbmc_scrublet_threshold_test.h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 1 : "output channel should contain one event, got ${output_list.size()}"
    assert output_list[0][0] == "pbmc_scrublet_threshold" : "Output id should be `pbmc_scrublet_threshold`."
  }

  test_ch = input_ch
    | workflow3_test.run(
      fromState: [
        "input": "output",
        "orig_input": "orig_input",
        "scrublet_score_threshold": "scrublet_score_threshold"
      ],
    )
}

workflow test_wf4 {

  resources_test = file(params.resources_test)

  // Shared single-sample processing parameters. QC metrics are calculated
  def base_args = [
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
    var_name_ribosomal_genes: 'ribosomal',
    obs_name_ribosomal_fraction: 'fraction_ribosomal',
    add_id_to_obs: true,
    add_id_make_observation_keys_unique: true,
    add_id_obs_output: "sample_id",
    intersect_obs: true
  ]

  input_ch = Channel.fromList([
    // Default output slots: QC metric columns use the default names.
    base_args + [
      id: "pbmc_default_slots",
      output: "pbmc_default_slots.h5mu"
    ],
    // Non-default output slots: QC metric columns use custom names and a
    // non-default set of top_n_vars. The renamed-away default names must not
    // appear in the output (verified by the test component).
    base_args + [
      id: "pbmc_custom_slots",
      output: "pbmc_custom_slots.h5mu",
      output_obs_num_nonzero_vars: "custom_num_nonzero_vars",
      output_obs_total_counts_vars: "custom_obs_total_counts",
      output_var_num_nonzero_obs: "custom_num_nonzero_obs",
      output_var_total_counts_obs: "custom_var_total_counts",
      output_var_obs_mean: "custom_obs_mean",
      output_var_pct_dropout: "custom_pct_dropout",
      top_n_vars: [10, 20]
    ]
  ])
  | map{ state -> [state.id, state] }
  | process_singlesample.run(
    toState: { id, output, state -> state + output + [orig_input: state.input] }
  )

  assert_ch = input_ch
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 2 : "output channel should contain two events, got ${output_list.size()}"
    assert output_list.collect{it[0]}.sort() == ["pbmc_custom_slots", "pbmc_default_slots"] : "Output ids should be `pbmc_custom_slots` and `pbmc_default_slots`."
  }

  test_ch = input_ch
    | workflow4_test.run(
      // Only pass slot arguments that were explicitly set on the event. Passing
      // a null would override the test component's default (which mirrors the
      // process_singlesample default), so unset arguments are omitted instead.
      fromState: { id, state ->
        def args = [
          "input": state.output,
          "orig_input": state.orig_input
        ]
        def passthrough = [
          "output_obs_num_nonzero_vars",
          "output_obs_total_counts_vars",
          "output_var_num_nonzero_obs",
          "output_var_total_counts_obs",
          "output_var_obs_mean",
          "output_var_pct_dropout",
          "log1p_transform",
          "top_n_vars",
          "var_name_mitochondrial_genes",
          "var_name_ribosomal_genes"
        ]
        passthrough.each { key ->
          if (state[key] != null) { args[key] = state[key] }
        }
        return args
      }
    )
}
