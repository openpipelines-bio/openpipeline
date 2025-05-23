nextflow.enable.dsl=2

include { cellranger_postprocessing } from params.rootDir + "/target/_test/nextflow/workflows/ingestion/cellranger_postprocessing/main.nf"
include { from_10xh5_to_h5mu } from params.rootDir + "/target/nextflow/convert/from_10xh5_to_h5mu/main.nf"
include { cellranger_postprocessing_test } from params.rootDir + "/target/nextflow/test_workflows/ingestion/cellranger_postprocessing_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "foo",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        input_og: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        perform_correction: true,
        min_genes: 100,
        min_counts: 1000,
        cellbender_epochs: 5
      ] 
    ])
    | map{ state -> [state.id, state] }

    | cellranger_postprocessing.run(
      toState: {id, output, state ->
        output + [
          input_og: state.input_og,
          perform_correction: state.perform_correction
        ]
      }
    )

    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }

    | cellranger_postprocessing_test.run(
      fromState: {id, state ->
        [
          input: state.output,
          input_og: state.input_og,
          is_corrected: state.perform_correction
        ]
      }
    )

    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}

workflow test_wf2 {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "zing",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        input_og: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        perform_correction: false,
        min_genes: 100,
        min_counts: 1000,
        cellbender_epochs: 5
      ]
    ])
    | map{ state -> [state.id, state] }
    | cellranger_postprocessing.run(
      toState: {id, output, state ->
        output + [
          input_og: state.input_og,
          perform_correction: state.perform_correction
        ]
      }
    )
    
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, out]"
      assert output[1] instanceof Map : "Output should be a Map."
      // todo: check whether output dir contains fastq files
      "Output: $output"
    }

    | cellranger_postprocessing_test.run(
      fromState: {id, state ->
        [
          input: state.output,
          input_og: state.input_og,
          is_corrected: state.perform_correction
        ]
      }
    )

    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
    }
}
