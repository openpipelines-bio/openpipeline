nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { add_id } from targetDir + "/metadata/add_id/main.nf"
include { split_modalities } from targetDir + '/dataflow/split_modalities/main.nf'
include { merge } from targetDir + '/dataflow/merge/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'
include { run_wf as rna_singlesample } from workflowDir + '/multiomics/rna_singlesample/main.nf'
include { run_wf as rna_multisample } from workflowDir + '/multiomics/rna_multisample/main.nf'
include { run_wf as integration } from workflowDir + '/multiomics/integration/main.nf'

include { readConfig; viashChannel; helpMessage; readCsv } from workflowDir + "/utils/WorkflowHelper.nf"
include { passthroughMap as pmap; passthroughFlatMap as pFlatMap } from workflowDir + "/utils/DataflowHelper.nf"
config = readConfig("$projectDir/config.vsh.yaml")


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
    // add ids to obs_names and to .obs[sample_id]
    | map { id, data ->
      def new_data = [
        input: data.input, 
        input_id: id, 
        make_observation_keys_unique: true, 
        obs_output: 'sample_id'
      ]
      def new_passthrough = ["obs_covariates": data.obs_covariates]
      [ id, new_data, new_passthrough ]
    }
    | add_id 

    // split by modality
    | split_modalities
    
    // combine output types csv
    | pFlatMap {id, data, passthrough ->
      def outputDir = data.output
      def types = readCsv(data.output_types.toString())
      
      types.collect{ dat ->
        // def new_id = id + "_" + dat.name
        def new_id = id // it's okay because the channel will get split up anyways
        def new_data = outputDir.resolve(dat.filename)
        def new_passthrough = passthrough + [ modality: dat.name ]
        [ new_id, new_data, new_passthrough]
      }
    }

    modality_processors = [
      [id: "rna", singlesample: rna_singlesample, multisample: rna_multisample],
      [id: "vdj_t", singlesample: null, multisample: null],
      [id: "vdj_b", singlesample: null, multisample: null],
      [id: "prot", singlesample: null, multisample: null],
    ]

    mod_chs = modality_processors.collect{ modality_processor ->
      mod_ch = start_ch
        | filter{ it[2].modality == modality_processor.id }

      ss_ch = (modality_processor.singlesample ? mod_ch | modality_processor.singlesample : mod_ch)

      input_ms_ch = ss_ch
        | toSortedList{ a, b -> b[0] <=> a[0] }
        | filter { it.size() != 0 } // filter when event is empty
        | map{ list -> 
          new_data = ["sample_id": list.collect{it[0]}, "input": list.collect{it[1]}]
          ["combined_" + modality_processor.id, new_data] + list[0].drop(2)
        }
      return (modality_processor.multisample ? \
              input_ms_ch | modality_processor.multisample : \
              input_ms_ch | concat.run(key: "concat_" + modality_processor.id,
                                       renameKeys: [input_id: "id"]))
    }


  output_ch = mod_chs[0].concat(*mod_chs.drop(1))
    | toSortedList{ a, b -> b[0] <=> a[0] }
    | map { list -> 
      ["merged", list.collect{it[1]}] + list[0].drop(2)
    }
    | merge
    | map { id, data, passthrough -> 
      new_data = [
        "input": data, 
        "layer": "log_normalized", 
        "obs_covariates": passthrough.get("obs_covariates")
      ]
      [id, new_data]
    }
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
    param_list: [
        [
          id: "mouse",
          input: params.resources_test + "/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
          obs_covariates: "sample_id",
          publish_dir: "foo/"
        ],
        [
          id: "human",
          input: params.resources_test + "/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
          input_type: "10xmtx",
          obs_covariates: "sample_id",
          publish_dir: "foo/"
        ]
      ]
    ]


  output_ch =
    viashChannel(testParams, config)
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
