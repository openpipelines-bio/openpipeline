nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { bd_rhapsody_wta } from targetDir + "/mapping/bd_rhapsody_wta/main.nf"
include { from_bdrhap_to_h5mu } from targetDir + "/convert/from_bdrhap_to_h5mu/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/ingestion/bd_rhapsody_wta/config.vsh.yaml")

// keep track of whether this is an integration test or not
global_params = [ do_publish: true ]

workflow {
  helpMessage(config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // Step 1: group fastq files per lane
    | flatMap { tup ->
      id = tup[0]
      data = tup[1].clone()

      // preproc input
      input = data.remove("input")
      if (input instanceof Path) {
        input = [ input ]
      }

      // get read regex
      r1r2_regex = data.remove("r1r2_regex")

      input_with_new_ids = input.collect { file ->
        new_id = file.name.replaceAll(r1r2_regex, '.').replaceAll("\\.fastq\\.gz", "")
        [ new_id, file ]
      }
      new_ids = input_with_new_ids.collect{it[0]}.unique()
      new_ids.collect { new_id -> 
        new_input = input_with_new_ids.findAll{it[0] == new_id}.collect{it[1]}
        assert new_input.size() == 2 : "Number of fastqs for id '$new_id' should be two. Found: ${new_input}"

        [ new_id, [ input: new_input ] + data, [ tuple_orig_id: id ] ]
      }
    }

    // Step 2: run BD rhapsody WTA
    | view { "running_bd_rhapsody: $it (orig_id: ${it[2].tuple_orig_id})" }
    | bd_rhapsody_wta.run(
      auto: [ publish: global_params.do_publish ]
    )

    // Step 3: group outputs per sample
    | map { id, input, extra -> [ extra.tuple_orig_id, input ] }
    | groupTuple()
    | map { id, inputs -> [ id, [ id: id, input: inputs] ] }

    // Step 4: convert to h5ad
    | view { "converting_to_h5mu: $it" }
    | from_bdrhap_to_h5mu.run(
      auto: [ publish: global_params.do_publish ]
    )

  emit:
  output_ch
}

workflow test_wf {
  // don't publish output
  global_params.do_publish = false

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    input: params.resources_test + "/bdrhap_5kjrt/raw/12WTA*.fastq.gz",
    reference_genome: params.resources_test + "/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz",
    transcriptome_annotation: params.resources_test + "/bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf",
    override_min_cores: 1,
    override_min_ram: 2,
    putative_cell_call: "mRNA",
    exact_cell_count: 4900
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
      | toList()
      | map { output_list ->
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
      }
      // | check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}