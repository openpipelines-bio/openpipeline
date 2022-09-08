nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { bd_rhapsody } from targetDir + "/mapping/bd_rhapsody/main.nf"
include { from_bdrhap_to_h5mu } from targetDir + "/convert/from_bdrhap_to_h5mu/main.nf"

include { readConfig; viashChannel; helpMessage } from workflowDir + "/utils/WorkflowHelper.nf"

config = readConfig("$workflowDir/ingestion/bd_rhapsody/config.vsh.yaml")

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

    | bd_rhapsody.run(auto: [ publish: true ])

    | map { id, file -> [ id, [ input: file, id: id ]]}
    | from_bdrhap_to_h5mu.run(auto: [ publish: true ])

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    id: "foo",
    mode: "wta",
    input: params.resources_test + "/bdrhap_5kjrt/raw/12WTA*.fastq.gz",
    reference: params.resources_test + "/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz",
    transcriptome_annotation: params.resources_test + "/bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf",
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
      | view { output_list ->
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "foo" : "Output ID should be same as input ID"
      }
      // | check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}