workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Avoid conflict with other output arguments
    | map {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | totalvi.run(
      fromState: {id, state -> [
        "input": state.input,
        "rna_modality": state.rna_modality,
        "prot_modality": state.prot_modality,
        "input_layer_rna": state.input_layer_rna,
        "input_layer_protein": state.input_layer_protein,
        "obs_batch": state.obs_batch,
        "obs_size_factor": state.obs_size_factor,
        "obs_categorical_covariate": state.obs_categorical_covariate,
        "obs_continuous_covariate": state.obs_continuous_covariate,
        "var_gene_names": state.var_gene_names,
        "var_protein_names": state.var_protein_names,
        "var_input": state.var_input,
        "max_epochs": state.max_epochs,
        "early_stopping": state.early_stopping,
        "obsm_output": state.obsm_integrated,
        "obsm_normalized_rna_output": state.obsm_normalized_rna_output,
        "obsm_normalized_protein_output": state.obsm_normalized_protein_output,
        "output_model": state.output_model
      ]},
      toState: [
        "input": "output",
        "output_model": "output_model"
      ]
    )

    | neighbors_leiden_umap.run(
      fromState: {id, state -> [
        "input": state.input,
        "modality": state.rna_modality,
        "uns_neighbors": state.uns_neighbors,
        "obsp_neighbor_distances": state.obsp_neighbor_distances,
        "obsp_neighbor_connectivities": state.obsp_neighbor_connectivities,
        "obs_cluster": state.obs_cluster,
        "leiden_resolution": state.leiden_resolution,
        "obsm_umap": state.obsm_umap,
        "obsm_input": state.obsm_integrated,
        "output": state.workflow_output
      ]},
      toState: ["output": "output"],
    )

    | setState(["output", "output_model"])
  emit:
  output_ch
}
