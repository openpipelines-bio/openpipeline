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
          aggregation_method: "aggregation_method",
          random_state: "random_state"
        ],
        args:[
          obs_cell_count: "n_cells"
        ],
        toState: [ "input": "output" ]
      )

      | delimit_counts.run(
          fromState: [
            input: "input",
            modality: "modality",
            min_counts: "min_obs_per_sample",

          ],
          args: [
            do_subset: "True",
            obs_count_column: "n_cells",
            obs_name_filter: "delimit_samples_per_pseudobulk",
          ],
          toState: [ "input": "output" ]
      )

      | filter_with_pattern.run(
          runIf: { id, state ->
            state.filter_gene_patterns
          },
          fromState: [
            input: "input",
            modality: "modality",
            var_gene_names: "var_gene_names",
            pattern: "filter_gene_patterns",
            min_counts: "min_obs_per_sample",

          ],
          args: [
            do_subset: "True",
            var_name_filter: "filter_with_gene_pattern"
          ],
          toState: [ "input": "output" ]
      )

      | filter_with_counts.run(
          fromState: [
            input: "input",
            modality: "modality",
            layer: "input_layer",
            min_cells_per_gene: "filter_genes_min_samples",
          ],
          args: [
            do_subset: "True",
            var_name_filter: "filter_with_counts"
          ],
          toState: [ "input": "output" ]
      )

      | deseq2.run(
        fromState: [
          input: "input",
          modality: "modality",
          input_layer: "input_layer",
          var_gene_names: "var_gene_names",
          obs_cell_group: "obs_cell_group",
          design_formula: "design_formula",
          contrast_column: "contrast_column",
          contrast_values: "contrast_values",
          p_adj_threshold: "p_adj_threshold",
          log2fc_threshold: "log2fc_threshold",
          output_dir: "output",
          output_prefix: "output_prefix"
        ],
        toState: [ "output": "output_dir" ]
      )

      | setState([ "output" ])

  emit:
    output_ch
}