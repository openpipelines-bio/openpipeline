nextflow.enable.dsl=2

workflowDir = params.rootDir + "/src/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { readConfig; helpMessage; channelFromParams; preprocessInputs } from workflowDir + "/utils/WorkflowHelper.nf"
include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"
include { do_filter } from targetDir + "/filter/do_filter/main.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap} from workflowDir + "/utils/DataflowHelper.nf"
include { qc as unfiltered_counts_qc_metrics_prot } from workflowDir + "/qc/qc/main.nf"

config = readConfig("$workflowDir/multiomics/prot_singlesample/config.vsh.yaml")

workflow prot_singlesample {
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
  output_ch = input_ch
    // split params for downstream components
    | preprocessInputs("config": config)
    | pmap {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
          
    | pmap { id, orig_state ->
      def new_state = [
        "input": orig_state.input,
        "var_name_mitondrial": null,
        "mitochondrial_gene_regex": null,
        "top_n_vars": orig_state.top_n_vars,
        "output": orig_state.output,
        "modality": "prot",
        "layer": null,
        "var_qc_metrics": orig_state.var_qc_metrics
      ]
      [id, new_state, orig_state]
    }

    | unfiltered_counts_qc_metrics_prot
    | pmap { id, state, orig_state ->
      stateMapping = [
          "min_counts": "min_counts",
          "max_counts": "max_counts",
          "min_genes_per_cell": "min_proteins_per_cell",
          "max_genes_per_cell": "max_proteins_per_cell",
          "min_cells_per_gene": "min_cells_per_protein",
          "workflow_output": "workflow_output"
      ]
      def new_state = stateMapping.collectEntries{newKey, origKey ->
        [newKey, orig_state[origKey]]
      }
      [id, new_state + ["input": state.output]]
    }

    // filtering
    | filter_with_counts.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "min_counts": state.min_counts,
          "max_counts": state.max_counts,
          "min_genes_per_cell": state.min_proteins_per_cell,
          "max_genes_per_cell": state.max_proteins_per_cell,
          "min_cells_per_gene": state.min_cells_per_protein,
          "obs_name_filter": "filter_with_counts",
          "var_name_filter": "filter_with_counts",
          "modality": "prot"
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      fromState : { id, state ->
        [
          "input": state.input,
          "obs_filter": "filter_with_counts",
          "modality": "prot",
          "var_filter": "filter_with_counts",
          "output_compression": "gzip",
          "output": state.workflow_output
        ]
      },
      toState: ["output": "output"],
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
    max_counts: 100000,
    min_genes_per_cell: 2,
    max_genes_per_cell: 10000,
    min_cells_per_gene: 10,
    min_fraction_mito: 0.2,
    max_fraction_mito: 0.8,
    output: "foo.final.h5mu",
  ]

  output_ch =
    channelFromParams(testParams, config)
    // Add test passthrough 
    | map {list -> list + [test_passthrough: "test"]}
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain two elements; [id, file, passthrough]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
      "Output: $output"
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
      assert (output_list.collect({it[1].output.getFileName().toString()}) as Set).equals(["foo.final.h5mu"] as Set)
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}