workflow run_wf {
  take:
  input_ch

  main:

  // Maps each modality to the sub-workflow arguments it needs, keyed by the sub-workflow
  // argument name and pointing to the top-level state key holding the value. Adding a new
  // modality is a matter of adding a qc_filter_<modality> sub-workflow and an entry here.
  def modality_arguments = [
    "rna": [
      "layer": "rna_layer",
      "min_count": "rna_min_count",
      "max_quantile": "rna_max_quantile",
      "mitochondrial_gene_regex": "mitochondrial_gene_regex",
      "var_gene_names": "var_gene_names",
      "skip_scrublet": "skip_scrublet",
      "scrublet_score_threshold": "scrublet_score_threshold",
      "scrublet_expected_doublet_rate": "scrublet_expected_doublet_rate",
      "scrublet_min_counts": "scrublet_min_counts",
      "scrublet_min_cells": "scrublet_min_cells",
      "scrublet_min_gene_variability_percent": "scrublet_min_gene_variability_percent",
      "scrublet_num_pca_components": "scrublet_num_pca_components",
    ],
    "prot": [
      "layer": "prot_layer",
      "min_count": "prot_min_count",
      "max_quantile": "prot_max_quantile",
    ],
  ].asImmutable()

  split_ch = input_ch
    // Set aside the workflow output target to avoid conflicts with component outputs.
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // Split the multimodal input into one file per modality for separate processing.
    | split_modalities.run(
      fromState: {id, state -> ["input": state.input]},
      toState: ["split_output": "output", "split_types": "output_types"]
    )
    // Expand the split directory + types CSV into one event per modality.
    | flatMap {id, state ->
      def outputDir = state.split_output
      def csv = state.split_types.splitCsv(strip: true, sep: ",").findAll{!it[0].startsWith("#")}
      def header = csv.head()
      def types = csv.tail().collect { row -> [header, row].transpose().collectEntries() }
      types.collect{ dat ->
        [id, state + ["input": outputDir.resolve(dat.filename), "modality": dat.name, "_meta": ["join_id": id]]]
      }
    }
    | view {"After splitting modalities: $it"}

  // Known modalities (rna, prot) get their QC + doublet keep-flags computed per modality.
  known_ch = split_ch
    | runEach(
      components: [qc_filter_rna, qc_filter_prot],
      filter: {id, state, component -> "qc_filter_" + state.modality == component.config.name},
      fromState: {id, state, component ->
        def args = modality_arguments.get(state.modality).collectEntries{key_, value_ -> [key_, state[value_]]}
        if (state.modality == "rna") {
          // delimit_fraction expects a fraction in [0, 1], not a percentage.
          args += ["max_fraction_mito": state.max_pct_counts_mt != null ? state.max_pct_counts_mt / 100.0 : null]
        }
        return args + ["id": id, "input": state.input]
      },
      toState: ["input": "output"]
    )

  // Modalities without a dedicated sub-workflow pass through unchanged.
  unknown_ch = split_ch
    | filter {id, state -> !modality_arguments.containsKey(state.modality)}

  output_ch = known_ch.mix(unknown_ch)
    // Regroup the per-modality files of each sample.
    | map {id, state -> [state._meta.join_id, state]}
    | groupTuple(by: 0, sort: "hash")
    | map { id, states ->
        def new_input = states.collect{it.input}
        def modalities = states.collect{it.modality}.unique()
        def other_state_keys = states.inject([].toSet()){ current_keys, state ->
          current_keys + state.keySet()
        }.minus(["input", "modality", "_meta"])
        def new_state = other_state_keys.inject([:]){ old_state, argument_name ->
          def argument_values = states.collect{it.get(argument_name)}.unique()
          assert argument_values.size() == 1, "Arguments should be the same across modalities. \
                                               Argument name: $argument_name, values: $argument_values"
          old_state + [(argument_name): argument_values[0]]
        }
        [id, new_state + ["input": new_input, "modalities": modalities]]
    }
    // Merge the per-modality files back into a single multimodal object. All keep-flags are
    // preserved in each modality's .obs; no cells have been dropped yet.
    | merge.run(
      fromState: ["input": "input"],
      args: [output_compression: "gzip"],
      toState: ["input": "output"]
    )
    // Per-sample cell-count report, computed from the still-unsubset keep-flags.
    | cell_count_report.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "prot_modality": state.modalities.contains("prot") ? "prot" : null,
          "sample_id_column": state.sample_id_column,
          "rna_filter_columns": ["filter_counts_rna", "filter_quantile_rna", "filter_mito_rna"],
          "scrublet_filter_column": state.skip_scrublet ? null : "filter_scrublet",
          "prot_filter_columns": state.modalities.contains("prot") ? ["filter_counts_prot", "filter_quantile_prot"] : null,
          "output": state.cell_count_report,
        ]
      },
      toState: {id, output, state -> state + ["cell_count_report_output": output.output]}
    )
    // Apply the RNA cell keep-flags.
    | do_filter.run(
      key: "rna_filter_cells",
      runIf: {id, state -> state.modalities.contains("rna")},
      fromState: {id, state ->
        def obs_filter = ["filter_counts_rna", "filter_quantile_rna", "filter_mito_rna"]
        if (!state.skip_scrublet) {
          obs_filter += ["filter_scrublet"]
        }
        [
          "input": state.input,
          "modality": "rna",
          "obs_filter": obs_filter,
        ]
      },
      toState: ["input": "output"]
    )
    // Flag and remove RNA genes expressed in too few cells (evaluated after cell filtering).
    | filter_with_counts.run(
      key: "rna_gene_filter",
      runIf: {id, state -> state.modalities.contains("rna")},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "layer": state.rna_layer,
          "var_name_filter": "filter_genes_rna",
          "min_cells_per_gene": state.filter_genes_min_cells,
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      key: "rna_filter_genes",
      runIf: {id, state -> state.modalities.contains("rna")},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "var_filter": ["filter_genes_rna"],
        ]
      },
      toState: ["input": "output"]
    )
    // Apply the protein cell keep-flags.
    | do_filter.run(
      key: "prot_filter_cells",
      runIf: {id, state -> state.modalities.contains("prot")},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "obs_filter": ["filter_counts_prot", "filter_quantile_prot"],
        ]
      },
      toState: ["input": "output"]
    )
    // Keep only cells passing every filter in every modality.
    | intersect_obs.run(
      runIf: {id, state -> state.intersect_obs && state.modalities.size() > 1},
      fromState: {id, state ->
        [
          "input": state.input,
          "modalities": state.modalities,
          "output": state.workflow_output,
        ]
      },
      toState: ["input": "output"]
    )
    | setState(["output": "input", "cell_count_report": "cell_count_report_output"])

  emit:
  output_ch
}
