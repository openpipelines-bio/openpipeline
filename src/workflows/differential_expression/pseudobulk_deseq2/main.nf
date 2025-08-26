workflow run_wf {
  take:
    input_ch

  main:

    output_ch = input_ch
      | create_pseudobulk.run(
        fromState: [
          id: "id",
          input: "input",
          modality: "modality",
          input_layer: "input_layer",
          obs_label: "obs_cell_group",
          obs_groups: "obs_groups",
          aggregation_method: "aggregation_method",
          min_obs_per_sample: "min_obs_per_sample",
          random_state: "random_state",
          obs_cell_count: "obs_cell_count"
        ],
        toState: [ "input": "output" ]
      )

      | filter_with_counts.run(
          fromState: [
            id: "id",
            input: "input",
            modality: "modality",
            layer: "input_layer",
            min_counts: "min_obs_per_sample",

            obs_cell_count: "obs_cell_count"
          ],
          args: [
            do_subset: "True",

          ]
          toState: [ "input": "output" ]
      )

      | do_filter.run(

      )

      | filter_genes_by_pattern.run(

      )

      | do_filter.run(

      )

      | deseq2.run(
        fromState: [
          id: "id",
          input: "input",
          modality: "modality",
          input_layer: "input_layer",
          var_gene_name: "var_gene_name",
          obs_cell_group: "obs_cell_group",
          design_formula: "design_formula",
          contrast_column: "contrast_column",
          contrast_values: "contrast_values",
          filter_genes_min_samples: "filter_genes_min_samples",
          p_adj_threshold: "p_adj_threshold",
          log2fc_threshold: "log2fc_threshold",
          filter_gene_patterns: "filter_gene_patterns",
          output_dir: "output",
          output_prefix: "output_prefix"
        ],
        toState: [ "output": "output_dir" ]
      )

      | setState([ "output" ])

  emit:
    output_ch
}