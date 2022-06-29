nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { bd_rhapsody_wta } from targetDir + "/mapping/bd_rhapsody_wta/main.nf"
include { bd_rhapsody_wta_1_10_1 } from targetDir + "/mapping/bd_rhapsody_wta_1_10_1/main.nf"
include { from_bdrhap_to_h5mu } from targetDir + "/convert/from_bdrhap_to_h5mu/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/viash_workflow_helper.nf"

params.bd_rhapsody_version = "1.9.1"

workflow {
  params.testing = false

  helpMessage(params, config)

  viashChannel(params, config)
    | run_wf
}

workflow run_wf {
  take:
  input_ch

  main:
  def bd_rhap_wta = bd_rhapsody_wta
  if (params.bd_rhapsody_version == "1.10.1") {
    bd_rhap_wta = bd_rhapsody_wta_1_10_1
  }

  output_ch = input_ch
    // Step 1: group fastq files per lane
    | flatMap { tup ->
      id = tup[0]
      data = tup[1]

      // preproc input
      input = data.remove("input")
      if (input instanceof Path) {
        input = [ input ]
      }

      input_with_new_ids = input.collect { file ->
        new_id = file.name.replaceAll("[^a-zA-Z0-9]R[12]_*\\.fastq\\.gz\$", "")
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
    | bd_rhap_wta

    // Step 3: group outputs per sample
    | map { id, input, extra -> [ extra.tuple_orig_id, input ] }
    | groupTuple()
    | map { id, inputs -> [ id, [ id: id, input: inputs] ] }

    // Step 4: convert to h5ad
    | view { "converting_to_h5mu: $it" }
    | from_bdrhap_to_h5mu.run(
      auto: [ publish: ! params.testing ]
    )

  emit:
  output_ch
}

workflow test_wf {
  params.testing = true
  
  output_ch =
    Channel.value(
      [
        "foo",
        [
          input: file(params.rootDir + "/resources_test/bd_rhapsody_wta_test/raw/*.fastq.gz"),
          reference_genome: file(params.rootDir + "/resources_test/bd_rhapsody_wta_test/raw/GRCh38_primary_assembly_genome_chr1.tar.gz"),
          transcriptome_annotation: file(params.rootDir + "/resources_test/bd_rhapsody_wta_test/raw/gencode_v40_annotation_chr1.gtf"),
          override_min_cores: 1,
          override_min_ram: 2
        ]
      ]
    )
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain three elements; [id, file]"
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