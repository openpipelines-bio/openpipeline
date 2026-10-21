// beyond/subpopulation_clustering - annotated atlas -> atlas with subpopulation labels
//
// Input:  [ id, state ] with state.input an atlas h5mu carrying a broad cell-type label.
// Output: the same atlas with obs[--obs_label] added.
//
// Channel convention throughout: [ id, state_map ]

workflow run_wf {
  take:
  input_ch  // [ id, { input: atlas.h5mu, ..params.. } ]

  main:

  output_ch = input_ch

    // Preserve the final output filename so it survives each toState swap.
    | map { id, state ->
        [ id, state + [ "workflow_output": state.output ] ]
      }

    // -- 1. Split the atlas by broad cell type -------------------------------------
    | split_h5mu.run(
        args: [
          "output": "split_by_celltype",
          "output_files": "split_files.csv",
          "output_compression": "gzip",
        ],
        fromState: { id, state -> [
          "input":       state.input,
          "obs_feature": state.obs_cell_type,
        ]},
        toState: { id, output, state ->
          state + [
            "split_output_dir":   output.output,
            "split_output_files": output.output_files,
          ]
        }
      )

    // One event per cell type, read from the index CSV that split_h5mu writes.
    | flatMap { id, state ->
        def lines = state.split_output_files.readLines()
        // Under -stub the index file is empty; emit nothing rather than failing.
        if (lines.size() < 2) {
          return []
        }
        def header = lines[0].split(",")*.trim()
        lines.drop(1).findAll { !it.startsWith("#") && !it.isBlank() }.collect { line ->
          def values = line.split(",")*.trim()
          def entry  = [header, values].transpose().collectEntries()
          def subtype_id = id + "_" + entry.name
          [ subtype_id,
            state + [
              "input":     state.split_output_dir.resolve(entry.filename),
              "cell_type": entry.name,
              "_meta":     [ "join_id": state._meta?.join_id ?: id ],
            ]
          ]
        }
      }

    // -- 2. Per cell type: PCA -> neighbours -> Leiden ------------------------------
    | pca.run(
        key: "pca_celltype",
        args: [
          "obsm_output": "X_pca_celltype",
          "output_compression": "gzip",
          "overwrite": true,
        ],
        fromState: { id, state -> [
          "input":  state.input,
          "output": state.workflow_output,
        ]},
        toState: [ "input": "output" ]
      )

    | neighbors_leiden_umap.run(
        fromState: { id, state -> [
          "input":             state.input,
          "obsm_input":        "X_pca_celltype",
          "leiden_resolution": state.leiden_resolution,
          "obs_cluster":       state.obs_label,
          "output":            state.workflow_output,
        ]},
        args: [
          "uns_neighbors":                "neighbors",
          "obsp_neighbor_distances":      "distances",
          "obsp_neighbor_connectivities": "connectivities",
        ],
        toState: [ "input": "output" ]
      )

    // -- 3. Fan-in: concatenate the cell-type objects back into one atlas -----------
    | map { id, state -> [ state._meta.join_id, id, state ] }
    | groupTuple(by: 0, sort: "hash")
    | map { atlas_id, cell_type_ids, states ->
        [
          atlas_id,
          [
            "input":    states.collect { it.input },
            "input_id": cell_type_ids,
            "output":   states[0].workflow_output,
          ]
        ]
      }
    | concatenate_h5mu.run(
        key: "concatenate_celltypes",
        fromState: [
          "input":    "input",
          "input_id": "input_id",
          "output":   "output",
        ],
        toState: [ "output": "output" ]
      )
    | setState(["output"])

  emit:
  output_ch
}
