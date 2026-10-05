workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state ->
      [id, state + ["workflow_output": state.output]]
    }
    | geneformer_tokenize.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "layer",
        "var_gene_ids": "var_gene_ids",
        "sanitize_ensembl_ids": "sanitize_ensembl_ids",
        "gene_map": "gene_map",
        "var_gene_names": "var_gene_names",
        "allow_non_integer": "allow_non_integer",
        "model_version": "model_version",
        "obsm_output": "obsm_tokens",
      ],
      toState: ["input": "output"]
    )
    | geneformer_embeddings_extract.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_tokens",
        "model": "model",
        "model_version": "model_version",
        "emb_mode": "emb_mode",
        "emb_layer": "emb_layer",
        "forward_batch_size": "forward_batch_size",
        "obsm_output": "obsm_output",
      ],
      toState: ["input": "output"]
    )
    | neighbors_leiden_umap.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_output",
        "output": "workflow_output",
        "uns_neighbors": "uns_neighbors",
        "obsp_neighbor_distances": "obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "obsp_neighbor_connectivities",
        "leiden_resolution": "leiden_resolution",
        "obs_cluster": "obs_cluster",
        "obsm_umap": "obsm_umap",
      ],
      toState: ["output": "output"]
    )
    | setState(["output"])

  emit:
  output_ch
}
