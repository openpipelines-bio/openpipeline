nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { clr } from targetDir + '/transform/clr/main.nf'
include { add_id } from targetDir + '/metadata/add_id/main.nf'
include { calculate_qc_metrics } from targetDir + '/qc/calculate_qc_metrics/main.nf'
include { qc as prot_qc } from workflowDir + '/qc/qc/main.nf'
include { readConfig; helpMessage; channelFromParams; preprocessInputs } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/prot_multisample/config.vsh.yaml")

workflow prot_multisample {
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
  processed_input_ch = input_ch
    | preprocessInputs("config": config)

  processed_input_ch
    | toSortedList
    | map  { list ->
        found_output_files = list.collect{it[1].getOrDefault("output", null)}.unique()
        assert found_output_files.size() < 2, "The specified output file is not the same for all samples. Found: $found_output_files"
    }

  output_ch = processed_input_ch
    | clr.run(
      fromState: ["input": "input"],
      toState: ["input": "output"],
      args: [ 
        output_layer: "clr", 
        modality: "prot"
      ]
    )
    | pmap { id, state ->
      stateMapping = [
        "input": "input",
        "top_n_vars": "top_n_vars",
        "output": "output"
      ]
      def new_state = stateMapping.collectEntries{newKey, origKey ->
        [newKey, state[origKey]]
      }
      static_args = [
        "var_qc_metrics": null,
        "input_layer": null, // layer: null to use .X and not log transformed
        "modality": "prot",
        "var_name_mitochondrial_genes": null,
        "var_qc_metrics": null
      ]
      new_state = new_state + static_args
      [id, new_state]
    }
    | prot_qc


  emit:
  output_ch
}


workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "adt_samples",
    sample_id: "pbmc",
    input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
  ]

  output_ch =
    channelFromParams(testParams, config)
    | map {list -> list + [test_passthrough: "test"]}
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 3 : "outputs should contain three elements; [id, file, passthrough]"
      assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "adt_samples" : "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}
