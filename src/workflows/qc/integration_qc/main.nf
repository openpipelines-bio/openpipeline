workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Set aside the output for this workflow to avoid conflicts
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output_report]
      [id, new_state]
    }
    // Calculate integration qc metrics and add them to the MuData file
    | integration_metrics.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": state.modality,
          "obsm_output": state.obsm_output,
          "obsm_embeddings": state.obsm_embeddings,
          "obs_batch_label": state.obs_batch_label,
          "obs_cell_label": state.obs_cell_label,
          "obs_cluster": state.obs_cluster,
          "uns_neighbors": state.uns_neighbors,
          "obsp_neighbor_connectivities": state.obsp_neighbor_connectivities,
          "output": state.output,
          "output_compression": state.output_compression
       ]
      },

      toState: [
        "input": "output", 
      ],
    )
    // Generate visualizations and a report based on the metrics
    | integration_report.run(
      fromState: { id, state ->
      [
        "input": state.input,
        "modality": state.modality,
        "var_gene_names": state.var_gene_names,
        "obs_batch_label": state.obs_batch_label,
        "obs_cell_label": state.obs_cell_label,
        "uns_neighbors": state.uns_neighbors, 
        "obsm_umap": state.obsm_umap,
        "output_report": state.output_report,
        "output_umap_batch": state.output_umap_batch,
        "output_umap_label": state.output_umap_label,
        "output_metrics": state.output_metrics
      ]
      },
      auto: [ publish: true ],
      toState: { id, output, state ->
        [ 
          output_report: output.output_report,
          output_umap_batch: output.output_umap_batch,
          output_umap_label: output.output_umap_label,
          output_metrics: output.output_metrics
        ]
      }
    )

  emit:
  output_ch
}
