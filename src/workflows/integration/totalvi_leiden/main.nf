workflow neighbors_leiden_umap {
  take:
  integrated_ch

  main:
  neighbors_ch = integrated_ch
    | find_neighbors.run(
      fromState: [
        "input": "input",
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_input": "obsm_output", // use output from scvi as input for neighbors,
        "query_modality": "modality"
      ],
      toState: ["input": "output"]
    )

  with_leiden_ch = neighbors_ch
    | filter{list -> list[1].leiden_resolution}
    | leiden.run(
      fromState: [
        "input": "input",
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obsm_name": "obs_cluster",
        "resolution": "leiden_resolution",
        "query_modality": "modality",
      ],
      toState: ["input": "output"]
    )
    | move_obsm_to_obs.run(
      fromState: [
        "input": "input",
        "obsm_key": "obs_cluster",
        "query_modality": "modality",
      ],
      toState: ["input": "output"]
    )

  without_leiden_ch = neighbors_ch
    | filter{list -> !list[1].leiden_resolution}

  output_ch = with_leiden_ch.mix(without_leiden_ch)
    | umap.run(
      fromState: [
          "input": "input",
          "uns_neighbors": "uns_neighbors",
          "obsm_output": "obsm_umap",
          "query_modality": "modality",
        ],
      toState: ["output": "output"]
    )

  emit:
  output_ch
}

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
        "layer": "layer",
        "obs_batch": "obs_batch",
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
    | map { id, state -> // for gene expression
      stateMapping = [
        "input": "input",
        "uns_neighbors": "rna_uns_neighbors",
        "obsp_neighbor_distances": "rna_obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "rna_obsp_neighbor_connectivities",
        "obsm_output": "rna_obsm_output",
        "obs_cluster": "rna_obs_cluster",
        "leiden_resolution": "rna_leiden_resolution",
        "uns_neighbors": "rna_uns_neighbors",
        "obsm_umap": "obsm_umap",
        "modality": "modality"
      ]
      def new_state = stateMapping.collectEntries{newKey, origKey ->
        [newKey, state[origKey]]
      }
      [id, new_state, state]
    }
    | neighbors_leiden_umap
    | map { id, state, orig_state -> // for ADT
      stateMapping = [
        "uns_neighbors": "prot_uns_neighbors",
        "obsp_neighbor_distances": "prot_obsp_neighbor_distances",
        "obsp_neighbor_connectivities": "prot_obsp_neighbor_connectivities",
        "obsm_output": "prot_obsm_output",
        "obs_cluster": "prot_obs_cluster",
        "leiden_resolution": "prot_leiden_resolution",
        "uns_neighbors": "prot_uns_neighbors",
        "obsm_umap": "obsm_umap",
        "modality": "prot_modality",
        "workflow_output": "workflow_output",
        "query_model_path": "query_model_path",
        "reference_model_path": "reference_model_path"
      ]
      def new_state = stateMapping.collectEntries{newKey, origKey ->
        [newKey, orig_state[origKey]]
      }
      [id, new_state + ["input": state.output]]
    }
    | neighbors_leiden_umap
    | publish.run(
      fromState: { id, state -> [
          "input": state.output,
          "output": state.workflow_output,
          "compression": "gzip"
        ]
      },
      toState: ["output", "output"]
    )
    | setState(["output", "reference_model_path", "query_model_path"])
  emit:
  output_ch
}
