workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | cellphonedb_component.run(
        fromState: [
          "input": "input",
          "modality": "modality",
          "layer": "layer",
          "groupby": "groupby",
          "gene_symbol": "gene_symbol",
          "counts_data": "counts_data",
          "method": "method",
          "iterations": "iterations",
          "threshold": "threshold",
          "pvalue": "pvalue",
          "degs_file": "degs_file",
          "microenvs_file": "microenvs_file",
          "active_tfs_file": "active_tfs_file",
          "score_interactions": "score_interactions",
          "subsampling": "subsampling",
          "subsampling_log": "subsampling_log",
          "subsampling_num_pc": "subsampling_num_pc",
          "subsampling_num_cells": "subsampling_num_cells",
          "output": "output",
          "output_compression": "output_compression",
          "output_uns_means": "output_uns_means",
          "output_uns_deconvoluted": "output_uns_deconvoluted",
          "output_uns_deconvoluted_percents": "output_uns_deconvoluted_percents",
          "output_uns_pvalues": "output_uns_pvalues",
          "output_uns_significant_means": "output_uns_significant_means",
          "output_uns_relevant_interactions": "output_uns_relevant_interactions",
          "output_uns_interaction_scores": "output_uns_interaction_scores",
          "output_uns_cellsign_active_interactions": "output_uns_cellsign_active_interactions",
          "output_uns_cellsign_active_interactions_deconvoluted": "output_uns_cellsign_active_interactions_deconvoluted"
        ],
        toState: ["output": "output"]
      )
      | setState(["output"])

  emit:
    output_ch
}
