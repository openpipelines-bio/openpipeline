nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { split_modalities } from targetDir + '/split/split_modalities/main.nf'
include { merge } from targetDir + '/integrate/merge/main.nf'
include { concat } from targetDir + '/integrate/concat/main.nf'
include { run_wf as rna_singlesample } from workflowDir + '/multiomics/rna_singlesample/main.nf'
include { run_wf as rna_multisample } from workflowDir + '/multiomics/rna_multisample/main.nf'
include { run_wf as integration } from workflowDir + '/multiomics/integration/main.nf'

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/multiomics/full_pipeline/config.vsh.yaml")


workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf

}

workflow run_wf {
  take:
  input_ch

  main:
  start_ch = input_ch
    | split_modalities

  rna_ch = start_ch
    | map { id, output_dir -> 
      files_list = output_dir.listFiles({ file -> file.name.endsWith('_rna.h5mu') && !file.isDirectory() })
      assert files_list.size() == 1
      [ id, [ input: files_list.first() ] ]
    }
    | rna_singlesample
    | toSortedList()
    | map { tups -> tups.transpose()}
    | map {id, files -> [id.join(','), ["id": id, "input": files]]}
    | rna_multisample


  atac_ch = start_ch
    | map { id, output_dir -> 
        [ id, 
          output_dir.listFiles({ file -> file.name.endsWith('_atac.h5mu') && !file.isDirectory() })
        ] }
    | map { id, files_list -> assert files_list.size() == 1; [id, files_list.first()] }
    | toSortedList()
    | map {list -> ["combined_samples_atac", 
                      ["input": list.collect{it[1]},
                       "sample_names":  list.collect{it[0]}]
                    ]}
    | concat // concat will be integrated into process_atac_multisample in the future
  
  output_ch = rna_ch.concat(atac_ch)
    | toSortedList()
    | map {list -> ["merged", list.collect{it[1]}]}
    | merge
    | map {id, path -> [id, ["input": path, "layer": "log_normalized"]]}
    | integration

  emit:
  output_ch
}

workflow test_wf {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "mouse;human",
    input: params.resources_test + "/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu;" + params.resources_test + "/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu",
    publish_dir: "foo/"
  ]

  output_ch =
    viashChannel(testParams, config)
      | flatMap {id, input_parameters -> [input_parameters.id, input_parameters.input].transpose()}
      | view { "Input: $it" }
      | run_wf
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1]}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
      }
}
