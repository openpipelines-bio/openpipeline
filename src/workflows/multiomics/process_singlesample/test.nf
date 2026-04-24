nextflow.enable.dsl=2
targetDir = params.rootDir + "/target"

include { process_singlesample } from targetDir + "/nextflow/workflows/multiomics/process_singlesample/main.nf"
include { workflow_test } from targetDir + "/_test/nextflow/test_workflows/multiomics/process_singlesample/workflow_test/main.nf"
include { workflow2_test } from targetDir + "/_test/nextflow/test_workflows/multiomics/process_singlesample/workflow2_test/main.nf"
include { workflow3_test } from targetDir + "/_test/nextflow/test_workflows/multiomics/process_singlesample/workflow3_test/main.nf"

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
      id: "pbmc",
      input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
      rna_min_counts: 2,
      rna_max_counts: 1000000,
      rna_min_genes_per_cell: 1,
      rna_max_genes_per_cell: 1000000,
      rna_min_percentile_counts: 0.05,
      rna_max_percentile_counts: 0.95,
      rna_min_cells_per_gene: 1,
      rna_min_fraction_mito: 0.0,
      rna_max_fraction_mito: 1.0,
      prot_min_counts: 3,
      prot_max_counts: 1000000,
      prot_min_percentile_counts: 0.05,
      prot_max_percentile_counts: 0.95,
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
    | workflow3_test.run(
      fromState: [
        "input": "output",
        "orig_input": "orig_input"
      ],
    )
}