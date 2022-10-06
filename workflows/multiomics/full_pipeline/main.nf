nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { split_modalities } from targetDir + '/dataflow/split_modalities/main.nf'
include { merge } from targetDir + '/dataflow/merge/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'
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
    // Store obs_covariates for later use
    | map { id, data ->
      new_data = data.clone()
      new_data.removeAll {k, v -> k == 'obs_covariates'}
      [ id, new_data, ["obs_covariates": data.obs_covariates] ]
    }
    | split_modalities

  rna_ch = start_ch
    | map { id, output_dir, other_params -> 
      files_list = output_dir.listFiles({ file -> file.name.endsWith('_rna.h5mu') && !file.isDirectory() })
      assert files_list.size() == 1
      [ id, [ input: files_list.first(), "id": id], other_params]
    }
    | rna_singlesample
    | toSortedList({ a, b -> b[0] <=> a[0] })
    | map { tups -> tups.transpose()}
    | map {id, files, other_params -> [id.join(','), ["id": id, "input": files]] + other_params.first()}
    | rna_multisample



  atac_ch = start_ch
    | map { id, output_dir, other_params -> 
        [ id, 
          output_dir.listFiles({ file -> file.name.endsWith('_atac.h5mu') && !file.isDirectory() }),
          other_params
        ]}
    | map { id, files_list, other_params -> assert files_list.size() == 1; [id, files_list.first(), other_params] }
    | toSortedList({ a, b -> b[0] <=> a[0] })
    | map {list -> ["combined_samples_atac", 
                      ["input": list.collect{it[1]},
                       "input_id":  list.collect{it[0]}],
                      list.collect{it[2]}.first()
                    ]}
    | concat // concat will be integrated into process_atac_multisample in the future

  output_ch = rna_ch.concat(atac_ch)
    | toSortedList()
    | map {list -> ["merged", list.collect{it[1]}] + list.collect{it[2]}.first()}
    | merge
    | map {tup -> [tup[0], ["input": tup[1], "layer": "log_normalized", "obs_covariates": tup[2].get("obs_covariates")]]}
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
    input: params.resources_test + "/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu;" + params.resources_test + "/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
    obs_covariates: ["sample_id"],
    publish_dir: "foo/",

  ]

  output_ch =
    viashChannel(testParams, config)
      | flatMap {id, input_parameters -> [input_parameters.id,
                                          input_parameters.input,
                                          Collections.nCopies(input_parameters.id.size(), 
                                                              input_parameters.obs_covariates)
                                         ].transpose()}
      | map {id,  input_file, obs_covariates -> [id, ["id": id, "input": input_file, "obs_covariates": obs_covariates]]}
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
