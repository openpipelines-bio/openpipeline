workflow run_wf {
  take:
  input_ch

  main:
  selected_ch = input_ch
    | map {id, state ->
      [id, state + ["workflow_output": state.output]]
    }
    // 1. disease / healthy / unused labels from the clusterings
    | label_cells.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obs_cluster": "obs_cluster",
        "obs_reference": "obs_reference",
        "disease_clusters": "disease_clusters",
        "disease_reference_clusters": "disease_reference_clusters",
        "healthy_clusters": "healthy_clusters",
        "healthy_reference_clusters": "healthy_reference_clusters",
      ],
      toState: ["input": "output"]
    )
    // 2. QC, depth filter, sampling and batches
    | sample_cells.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "input_layer": "layer",
        "var_gene_names": "var_gene_names",
        "skip_qc": "skip_qc",
        "nmads": "nmads",
        "nmads_mt": "nmads_mt",
        "pct_mt_max": "pct_mt_max",
        "min_genes": "min_genes",
        "n_disease_cells": "n_disease_cells",
        "seed": "sample_seed",
        "batch_size": "batch_size",
      ],
      toState: ["input": "output"]
    )
    // 3. healthy and disease centroids; this file is the h5mu output
    | compute_centroids.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_input": "obsm_embedding",
        "output": "workflow_output",
      ],
      args: [
        "obs_group": "perturbation_group",
        "obs_filter": "perturbation_selected",
        "groups": ["disease", "healthy"],
      ],
      toState: ["input": "output", "output": "output"]
    )

  // 4. tokenize the cells to perturb, uncropped, and split them per batch
  batches_ch = selected_ch
    | do_filter.run(
      fromState: [
        "input": "input",
        "modality": "modality",
      ],
      args: ["obs_filter": ["perturbation_perturb"]],
      toState: ["perturb_h5mu": "output"]
    )
    | geneformer_tokenize.run(
      fromState: [
        "input": "perturb_h5mu",
        "modality": "modality",
        "input_layer": "layer",
        "var_gene_ids": "var_gene_ids",
        "sanitize_ensembl_ids": "sanitize_ensembl_ids",
        "gene_map": "gene_map",
        "var_gene_names": "var_gene_names",
        "model_version": "model_version",
      ],
      args: ["obsm_output_uncropped": "geneformer_tokens_uncropped"],
      toState: ["perturb_h5mu": "output"]
    )
    | split_h5mu.run(
      fromState: [
        "input": "perturb_h5mu",
        "modality": "modality",
      ],
      args: ["obs_feature": "perturbation_batch"],
      toState: ["batch_dir": "output", "batch_files": "output_files"]
    )
    // one event per batch
    | flatMap {id, state ->
      state.batch_files.splitCsv(header: true)
        .sort { it.name }
        .collect { row ->
          [
            "${id}_${row.name}".toString(),
            state + [
              "run_id": id,
              "batch_h5mu": state.batch_dir.resolve(row.filename),
            ]
          ]
        }
    }

  // 5. per batch: knockouts, embedding and similarity shift
  shifts_ch = batches_ch
    | geneformer_virtual_cells.run(
      fromState: [
        "input": "batch_h5mu",
        "modality": "modality",
        "var_gene_names": "var_gene_names",
        "top_n_genes": "top_n_genes",
      ],
      toState: ["virtual_h5mu": "output"]
    )
    | geneformer_embeddings_extract.run(
      fromState: [
        "input": "virtual_h5mu",
        "modality": "modality",
        "model": "model",
        "model_version": "model_version",
        "emb_mode": "emb_mode",
        "emb_layer": "emb_layer",
        "forward_batch_size": "forward_batch_size",
        "obsm_output": "obsm_embedding",
      ],
      toState: ["virtual_h5mu": "output"]
    )
    | similarity_shift.run(
      fromState: [
        "input": "virtual_h5mu",
        "original": "batch_h5mu",
        "modality": "modality",
        "obsm_input": "obsm_embedding",
      ],
      toState: ["shift_csv": "output"]
    )

  // 6. gather the batches and rank the genes
  output_ch = shifts_ch
    | map {id, state -> [state.run_id, state]}
    | groupTuple(sort: { a, b -> a.batch_h5mu.name <=> b.batch_h5mu.name })
    | map {id, states ->
      [id, states[0] + ["output_shifts": states.collect { it.shift_csv }]]
    }
    | rank_genes.run(
      fromState: [
        "input": "output_shifts",
        "shift_column": "shift_column",
        "min_coverage": "min_coverage",
        "n_random": "n_random",
        "seed": "rank_seed",
        "pvalue_cutoff": "pvalue_cutoff",
        "top_n": "top_n",
        "output_ranked": "output_ranked",
        "output_top": "output_top",
      ],
      toState: ["output_ranked": "output_ranked", "output_top": "output_top"]
    )
    | setState(["output", "output_ranked", "output_top", "output_shifts"])

  emit:
  output_ch
}
