nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { filter_with_counts } from targetDir + "/filter/filter_with_counts/main.nf"
include { filter_with_scrublet } from targetDir + "/filter/filter_with_scrublet/main.nf"
include { do_filter } from targetDir + "/filter/do_filter/main.nf"

include { readConfig; helpMessage; channelFromParams; preprocessInputs } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/rna_singlesample/config.vsh.yaml")

workflow {
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
    | preprocessInputs("config": config)
    | setWorkflowArguments(
      filter_with_counts: [
          "min_counts": "min_counts",
          "max_counts": "max_counts",
          "min_genes_per_cell": "min_genes_per_cell",
          "max_genes_per_cell": "max_genes_per_cell",
          "min_cells_per_gene": "min_cells_per_gene",
          "min_fraction_mito": "min_fraction_mito",
          "max_fraction_mito": "max_fraction_mito",
      ],
      do_filter: [:],
      filter_with_scrublet: [:]
    )

    // cell filtering
    | getWorkflowArguments(key: "filter_with_counts")
    | filter_with_counts.run(
        args: [
          var_gene_names: "gene_symbol",
          obs_name_filter: "filter_with_counts", 
          var_name_filter: "filter_with_counts"
        ]
    )
    | getWorkflowArguments(key: "do_filter")
    | do_filter.run(
      args: [
        obs_filter: "filter_with_counts",
        var_filter: "filter_with_counts"
      ]
    )
    // doublet calling
    | getWorkflowArguments(key: "filter_with_scrublet")
    | view {"scrublet: $it"}
    | filter_with_scrublet.run(
      auto: [ publish: true, output_compression: "gzip" ]
    )
    | map {list -> [list[0], list[1]] + list.drop(3)}
    // TODO: ambient rna correction

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
    max_fraction_mito: 0.2
  ]

  output_ch =
    channelFromParams(testParams, config)
    // Add test passthrough 
    | map {list -> list + [test_passthrough: "test"]} 
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
    | view { output ->
      assert output.size() == 3 : "outputs should contain two elements; [id, file, passthrough]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}