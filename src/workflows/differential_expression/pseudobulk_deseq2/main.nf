workflow run_wf {
  take:
    input_ch

  main:

    output_ch = input_ch
      | create_pseudobulk.run(
        fromState: [
          input: "input",
          modality: "modality",
          input_layer: "input_layer",
          obs_label: "obs_cell_group",
          obs_groups: "obs_groups",
          aggregation_method: "aggregation_method"
          random_state: "random_state"
        ],
        args:[
          obs_cell_count: "n_cells"
        ]
        toState: [ "input": "output" ]
      )

      | delimit_counts.run(
          fromState: [
            input: "input",
            modality: "modality",
            obs_count_column: "input_layer",
            min_counts: "min_obs_per_sample",

          ],
          args: [
            do_subset: "True",
            obs_cell_count: "n_cells",
            obs_name_filter: "delimit_samples_per_pseudobulk",
          ]
          toState: [ "input": "output" ]
      )
      | filter_genes_by_pattern.run(
          
      )
      | filter_with_counts.run(
          fromState: [
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



      | do_filter.run(

      )

      | deseq2.run(
        fromState: [
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