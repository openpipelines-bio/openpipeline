nextflow.enable.dsl=2

workflowDir = params.rootDir + "/src/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"
include { filter_with_scrublet } from targetDir + "/filter/filter_with_scrublet/main.nf"
include { do_filter } from targetDir + "/filter/do_filter/main.nf"
include { qc as unfiltered_counts_qc_metrics_rna } from workflowDir + "/qc/qc/main.nf"
include { delimit_fraction } from targetDir + '/filter/delimit_fraction/main.nf'

include { readConfig; helpMessage; channelFromParams; preprocessInputs } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/rna_singlesample/config.vsh.yaml")

workflow rna_singlesample {
  helpMessage(config)

  channelFromParams(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }

}

workflow run_wf {
  take:
  input_ch

  main:
  qc_channel = input_ch
    | preprocessInputs("config": config)
    // Set aside the output for this workflow to avoid conflicts
    | pmap {id, state -> 
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // Check for correctness of mitochondrial gene detection arguments
    | pmap { id, state ->
      if (state.obs_name_mitochondrial_fraction && !state.var_name_mitochondrial_genes) {
        throw new RuntimeException("Using --obs_name_mitochondrial_fraction requires --var_name_mitochondrial_genes.")
      }
      if ((state.min_fraction_mito  || state.max_fraction_mito) && !state.obs_name_mitochondrial_fraction) {
        throw new RuntimeException("Enabling --min_fraction_mito or --max_fraction_mito requires --obs_name_mitochondrial_fraction.")
      }
      if (state.var_gene_names && !state.var_name_mitochondrial_genes) {
        System.err.println("Warning: --var_gene_names is set, but not --var_name_mitochondrial_genes. \
                           --var_gene_names is only required for mitochondrial gene detection and does \
                           nothing while not also setting --var_name_mitochondrial_genes")
      }
      if (state.mitochondrial_gene_regex && !state.var_name_mitochondrial_genes) {
        System.err.println("Warning: --mitochondrial_gene_regex is set, but not --var_name_mitochondrial_genes. \
                           --mitochondrial_gene_regex is only required for mitochondrial gene detection and does \
                           nothing while not also setting --var_name_mitochondrial_genes")
      }
      [id, state]
    }
    | pmap { id, orig_state ->
      // The rna singlesample processing allows detecting mitochondrial genes and filtering based
      // on the fraction of mitochondrial genes per cell
      // This behaviour is optional based on the presence of var_name_mitochondrial_genes
      // The behavior of other components must be tuned to this argument as well

      def new_state = [
        "input": orig_state.input,
        "top_n_vars": orig_state.top_n_vars,
        "output": orig_state.output,
        "modality": "rna",
        "layer": null
      ]
      if (orig_state.var_name_mitochondrial_genes){
        // Check if user has defined var columns to calculate metrics
        def new_var_qc_metrics = orig_state.var_qc_metrics != null ? orig_state.var_qc_metrics : []
        assert new_var_qc_metrics instanceof List
        // Add the mitochondrial genes var column to the columns to calculate statistics for if set.
        new_var_qc_metrics = ((new_var_qc_metrics as Set) + [orig_state.var_name_mitochondrial_genes]) as List

        def fraction_column_name = orig_state.obs_name_mitochondrial_fraction ? orig_state.obs_name_mitochondrial_fraction : "fraction_$orig_state.var_name_mitochondrial_genes";
        new_state += [
          "var_qc_metrics": new_var_qc_metrics,
          "obs_name_mitochondrial_fraction": fraction_column_name,
          "var_gene_names": orig_state.var_gene_names,
          "var_name_mitochondrial_genes": orig_state.var_name_mitochondrial_genes,
          "mitochondrial_gene_regex": orig_state.mitochondrial_gene_regex
        ]
      }
    
      return [id, new_state, orig_state]
    }
    | unfiltered_counts_qc_metrics_rna
    | pmap { id, state, orig_state ->
      stateMapping = [
        "min_counts": "min_counts",
        "max_counts": "max_counts",
        "min_genes_per_cell": "min_genes_per_cell",
        "max_genes_per_cell": "max_genes_per_cell",
        "min_cells_per_gene": "min_cells_per_gene",
        "var_name_mitochondrial_genes": "var_name_mitochondrial_genes",
        "workflow_output": "workflow_output"
      ]
      if (orig_state.min_fraction_mito) {
        stateMapping += ["min_fraction_mito": "min_fraction_mito"]
      }
      if (orig_state.max_fraction_mito) {
        stateMapping += ["max_fraction_mito": "max_fraction_mito"] 
      }
      def new_state = stateMapping.collectEntries{newKey, origKey ->
        [newKey, orig_state[origKey]]
      }
      [id, new_state + ["input": state.output]]
    }

  with_delimit_ch = qc_channel
    | filter {it[1].var_name_mitochondrial_genes}
    | view {"Before delimit fraction $it"}
    | delimit_fraction.run(
      fromState: {id, state -> 
        [
          "input": state.input,
          "obs_name_filter": "filter_mitochondrial",
          "min_fraction": state.min_fraction_mito,
          "max_fraction": state.max_fraction_mito,
          "obs_fraction_column": "fraction_$state.var_name_mitochondrial_genes"
        ]
      },
      toState: ["input": "output"]
    )
  without_delimit_ch = qc_channel
    | filter {!it[1].var_name_mitochondrial_genes}
  
  output_ch = with_delimit_ch.mix(without_delimit_ch)
    // cell filtering
    | filter_with_counts.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "obs_name_filter": "filter_with_counts",
          "var_name_filter": "filter_with_counts",
          "min_counts": state.min_counts,
          "max_counts": state.max_counts,
          "min_genes_per_cell": state.min_genes_per_cell,
          "max_genes_per_cell": state.max_genes_per_cell,
          "min_cells_per_gene": state.min_cells_per_gene,
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      fromState: {id, state ->
        def stateMapping = [
          input: state.input,
          var_filter: ["filter_with_counts"]
        ]
        def obs_filter = ["filter_with_counts"]
        if (state.var_name_mitochondrial_genes) {
          obs_filter += ["filter_mitochondrial"]
        }
        stateMapping += ["obs_filter": obs_filter]
        return stateMapping
      },
      toState: ["input": "output"]
    )
    // doublet calling
    | filter_with_scrublet.run(
      fromState: [
        "input": "input",
        "output": "workflow_output"
      ],
      args: [output_compression: "gzip"],
      auto: [ publish: true ]
    )

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    min_counts: 3,
    max_counts: 10000000,
    min_genes_per_cell: 2,
    max_genes_per_cell: 10000000,
    min_cells_per_gene: 2,
    min_fraction_mito: 0.05,
    max_fraction_mito: 0.2,
    var_gene_names: "gene_symbol",
    var_name_mitochondrial_genes: "mitochondrial",
    obs_name_mitochondrial_fraction: "fraction_mitochondrial",
    output: "foo.final.h5mu"
  ]

  output_ch =
    channelFromParams(testParams, config)
    // Add test passthrough 
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["foo.final.h5mu"] as Set)

    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}

workflow test_wf2 {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    min_counts: 3,
    max_counts: 10000000,
    min_genes_per_cell: 2,
    max_genes_per_cell: 10000000,
    min_cells_per_gene: 2,
    output: "foo.final.h5mu"
  ]

  output_ch =
    channelFromParams(testParams, config)
    // Add test passthrough 
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["foo.final.h5mu"] as Set)

    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}