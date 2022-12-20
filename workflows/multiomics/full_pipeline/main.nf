nextflow.enable.dsl=2
workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { add_id } from targetDir + "/metadata/add_id/main.nf"
include { split_modalities } from targetDir + '/dataflow/split_modalities/main.nf'
include { merge } from targetDir + '/dataflow/merge/main.nf'
include { concat } from targetDir + '/dataflow/concat/main.nf'
include { remove_modality }  from targetDir + '/filter/remove_modality/main.nf'
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
      def new_passthrough = ["obs_covariates": data.obs_covariates, 
                             "filter_with_hvg_var_output": data.filter_with_hvg_var_output]
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
      // [id: "vdj_t", singlesample: null, multisample: null],
      // [id: "vdj_b", singlesample: null, multisample: null],
      // [id: "prot", singlesample: null, multisample: null],
      // [id: "atac", singlesample: null, multisample: null],
    ]
    known_modalities = modality_processors.collect{it.id}

    unknown_channel = start_ch
      | filter { ! known_modalities.contains(it[2].modality)}
      | map { lst ->
          [lst[2].modality] + lst
       }
      | groupTuple(by: 0)
      // [modality, [sample_id1, sample_id2, ...], [h5mu1, h5mu2, ...], [passthrough, copy_passthrough, ...]]
      | map { grouped_lst ->
        modality_name = grouped_lst[0]
        ["combined_$modality_name", ["sample_id": grouped_lst[1], 
         "input": grouped_lst[2]], grouped_lst[3][0]] // passthrough is copied, just pick the first
      }
      | concat.run(renameKeys: [input_id: "sample_id"])

    mod_chs = modality_processors.collect{ modality_processor ->
      // Select the files corresponding to the currently selected modality
      mod_ch = start_ch
        | filter{ it[2].modality == modality_processor.id }
        | view {"start channel-$modality_processor.id: $it"}

      // Run the single-sample processing if defined
      ss_ch = (modality_processor.singlesample ? \
               mod_ch | modality_processor.singlesample : \
               mod_ch)
      
      // Reformat arguments to serve to the multisample processing
      input_ms_ch = ss_ch
        | view { "single-sample-input-$modality_processor.id: $it" }
        | toSortedList{ a, b -> b[0] <=> a[0] }
        | filter { it.size() != 0 } // filter when event is empty
        | map{ list -> 
          new_data = ["sample_id": list.collect{it[0]}, "input": list.collect{it[1]}]
          ["combined_$modality_processor.id", new_data] + list[0].drop(2)
        }
        | view { "input multichannel-$modality_processor.id: $it" }
      
      // Run the multisample processing if defined, otherwise just concatenate samples together
      out_ch = (
        modality_processor.multisample ? \
          input_ms_ch | modality_processor.multisample : \
          input_ms_ch | concat.run(
            key: "concat_" + modality_processor.id,
            renameKeys: [input_id: "sample_id"]
          )
      )
      return out_ch
    }

  // Concat channel if more than one modality was found
  known_mods_ch = mod_chs.size() > 1? mod_chs[0].concat(*mod_chs.drop(1)): mod_chs[0]
  merge_ch = unknown_channel.concat(known_mods_ch)
  output_ch = merge_ch
    | toSortedList{ a, b -> b[0] <=> a[0] }
    | map { list -> 
      ["merged", list.collect{it[1]}] + list[0].drop(2)
    }
    | merge
    | map { id, data, passthrough -> 
      new_data = [
        "input": data, 
        "layer": "log_normalized", 
        "obs_covariates": passthrough.get("obs_covariates"),
        "filter_with_hvg_var_output": passthrough.get("filter_with_hvg_var_output")
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
  input_ch = viashChannel(testParams, config)
    // store output value in 3rd slot for later use
    // and transform for concat component
    | map { tup ->
      data = tup[1]
      new_data = [ id: data.id, input: data.input ]
      [tup[0], new_data, data] + tup.drop(2)
    }

    human_ch = input_ch
      | filter{it[0] == "human"}
      | remove_modality.run(
        args: [ modality: "atac" ]
      )

    mouse_ch = input_ch
      | filter{it[0] == "mouse"}
      | remove_modality.run(
        args: [ modality: "rna" ]
      )

    output_ch_test_2 = human_ch.concat(mouse_ch)
      // Put back values for other parameters after removing modalities
      | map { tup ->
        tup[2].remove("input")
        new_data = [input: tup[1]]
        new_data.putAll(tup[2])
        [tup[0], new_data]
      }
      | view { "Input test 2: $it" }
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
