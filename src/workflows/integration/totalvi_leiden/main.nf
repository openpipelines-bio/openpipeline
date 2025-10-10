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
      fromState: [
        "input": "input",
        "input_layer": "layer",
        "obs_batch": "obs_batch",
        "obs_size_factor": "obs_size_factor",
        "obs_categorical_covariate": "obs_categorical_covariate",
        "obs_continuous_covariate": "obs_continuous_covariate",
        "query_modality": "modality",
        "query_proteins_modality": "prot_modality",
        "query_model_path": "query_model_path",
        "obsm_normalized_rna_output": "rna_obsm_output",
        "obsm_normalized_protein_output": "prot_obsm_output",
        "reference_model_path": "reference_model_path",
        "reference_modality": "rna_reference_modality",
        "reference_proteins_modality": "prot_reference_modality",
        "var_input": "var_input",
        "force_retrain": "force_retrain",
        "weight_decay": "weight_decay",
        "max_epochs": "max_epochs",
        "max_query_epochs": "max_query_epochs",
        "reference": "reference"
      ],
      toState: [
        "input": "output",
        "query_model_path": "query_model_path",
        "reference_model_path": "reference_model_path",
      ]
    )
    | neighbors_leiden_umap.run( // For gene expression
      key: "rna_neighbors_leiden_umap",
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "rna_obsm_output",
        "uns_neighbors": "rna_uns_neighbors",
        "obsp_neighbor_distances": "rna_obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "rna_obsp_neighbor_connectivities",
        "obs_cluster": "rna_obs_cluster",
        "leiden_resolution": "rna_leiden_resolution",
        "uns_neighbors": "rna_uns_neighbors",
        "obsm_umap": "obsm_umap",
      ],
      toState: ["input": "output"],
    )
    | neighbors_leiden_umap.run( // For ADT
      key: "adt_neighbors_leiden_umap",
      fromState: [
        "input": "input",
        "modality": "prot_modality",
        "obsm_input": "prot_obsm_output",
        "uns_neighbors": "prot_uns_neighbors",
        "obsp_neighbor_distances": "prot_obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "prot_obsp_neighbor_connectivities",
        "obs_cluster": "prot_obs_cluster",
        "leiden_resolution": "prot_leiden_resolution",
        "uns_neighbors": "prot_uns_neighbors",
        "obsm_umap": "obsm_umap",
        "output": "workflow_output",
      ],
      toState: ["output": "output"],
    )
    | setState(["output", "reference_model_path", "query_model_path"])
  emit:
  output_ch
}
