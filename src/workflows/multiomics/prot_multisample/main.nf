workflow run_wf {
  take:
  input_ch

  main:

  input_ch
    | toSortedList
    | map { list ->
        found_output_files = list.collect{it[1].getOrDefault("output", null)}.unique()
        assert found_output_files.size() < 2, "The specified output file is not the same for all samples. Found: $found_output_files"
    }

  output_ch = input_ch
    | clr.run(
      fromState: ["input": "input"],
      toState: ["input": "output"],
      args: [ 
        output_layer: "clr", 
        modality: "prot"
      ]
    )
    | prot_qc.run(
      key: "prot_qc",
      fromState: { id, state ->
        def newState = [
          "id": id,
          "input": state.input,
          "top_n_vars": state.top_n_vars,
          "var_qc_metrics": null,
          "input_layer": null, // layer: null to use .X and not log transformed
          "modality": "prot",
          "var_name_mitochondrial_genes": null,
          "var_qc_metrics": null
        ]
        newState
      }
    )
    | setState(["output"])


  emit:
  output_ch
}